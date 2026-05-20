#!/usr/bin/env nextflow

/*
 * SNPhylo2 Nextflow Pipeline
 * 
 * A scalable, reproducible pipeline for phylogenomic inference
 */

nextflow.enable.dsl=2

// Pipeline parameters
params.input = null
params.output = "results"
p params.config = null
params.threads = 4

// Validate inputs
if (!params.input) {
    error "Parameter 'input' is required. Specify with --input <vcf_file>"
}

// Parse input
input_path = file(params.input)
if (!input_path.exists()) {
    error "Input file does not exist: ${params.input}"
}

// Workflow definition
workflow SNPhylo2 {
    take:
        vcf_ch
    
    main:
        // Step 1: Quality Control
        QC(vcf_ch)
        
        // Step 2: Filtering
        FILTER(vcf_ch)
        
        // Step 3: LD Pruning
        PRUNE(FILTER.out.filtered_vcf)
        
        // Step 4: Tree Building
        TREE(PRUNE.out.pruned_vcf)
        
        // Step 5: Report Generation
        REPORT(
            QC.out.qc_report,
            FILTER.out.filter_stats,
            PRUNE.out.prune_stats,
            TREE.out.tree_file,
            TREE.out.tree_stats
        )
    
    emit:
        tree = TREE.out.tree_file
        report = REPORT.out.report
}

// Individual processes

process QC {
    tag "${vcf.name}"
    label 'low_mem'
    
    input:
        path vcf
    
    output:
        path "qc_report.html", emit: qc_report
        path "qc_stats.json", emit: qc_stats
    
    script:
        """
        snphylo2 qc \\
            -v ${vcf} \\
            -o qc_report.html
        
        # Extract stats for downstream use
        echo '{"step": "qc", "status": "complete"}' > qc_stats.json
        """
    
    stub:
        """
        touch qc_report.html
        echo '{"step": "qc", "status": "stub"}' > qc_stats.json
        """
}

process FILTER {
    tag "${vcf.name}"
    label 'medium_mem'
    
    input:
        path vcf
    
    output:
        path "*.filtered.vcf.gz", emit: filtered_vcf
        path "filter_stats.json", emit: filter_stats
    
    script:
        def prefix = vcf.name.replace('.vcf.gz', '').replace('.vcf', '')
        """
        snphylo2 filter \\
            -v ${vcf} \\
            -o ${prefix}.filtered.vcf.gz \\
            --maf 0.05 \\
            --max-missing 0.2 \\
            --min-depth 5
        
        # Generate stats
        echo '{"step": "filter", "status": "complete"}' > filter_stats.json
        """
    
    stub:
        def prefix = vcf.name.replace('.vcf.gz', '').replace('.vcf', '')
        """
        touch ${prefix}.filtered.vcf.gz
        echo '{"step": "filter", "status": "stub"}' > filter_stats.json
        """
}

process PRUNE {
    tag "${vcf.name}"
    label 'high_mem'
    
    input:
        path vcf
    
    output:
        path "*.pruned.vcf.gz", emit: pruned_vcf
        path "prune_stats.json", emit: prune_stats
    
    script:
        def prefix = vcf.name.replace('.filtered.vcf.gz', '')
        """
        snphylo2 prune \\
            -i ${vcf} \\
            -o ${prefix}.pruned.vcf.gz \\
            --window 50 \\
            --step 10 \\
            --r2 0.2 \\
            -t ${task.cpus}
        
        echo '{"step": "prune", "status": "complete"}' > prune_stats.json
        """
    
    stub:
        def prefix = vcf.name.replace('.filtered.vcf.gz', '')
        """
        touch ${prefix}.pruned.vcf.gz
        echo '{"step": "prune", "status": "stub"}' > prune_stats.json
        """
}

process TREE {
    tag "${vcf.name}"
    label 'high_cpu'
    
    input:
        path vcf
    
    output:
        path "*.tree.nwk", emit: tree_file
        path "tree_stats.json", emit: tree_stats
    
    script:
        def prefix = vcf.name.replace('.pruned.vcf.gz', '')
        """
        snphylo2 tree \\
            -i ${vcf} \\
            -o ${prefix}.tree.nwk \\
            --engine iqtree2 \\
            --model GTR+ASC \\
            --bootstrap 1000 \\
            -t ${task.cpus}
        
        echo '{"step": "tree", "status": "complete"}' > tree_stats.json
        """
    
    stub:
        def prefix = vcf.name.replace('.pruned.vcf.gz', '')
        """
        echo '((A,B),(C,D));' > ${prefix}.tree.nwk
        echo '{"step": "tree", "status": "stub"}' > tree_stats.json
        """
}

process REPORT {
    publishDir params.output, mode: 'copy'
    
    input:
        path qc_report
        path filter_stats
        path prune_stats
        path tree_file
        path tree_stats
    
    output:
        path "final_report.html", emit: report
        path "pipeline_stats.json", emit: stats
    
    script:
        """
        # Combine all stats
        cat ${filter_stats} ${prune_stats} ${tree_stats} > pipeline_stats.json
        
        # Generate final report
        snphylo2 report \\
            -d . \\
            -o final_report.html
        """
    
    stub:
        """
        touch final_report.html
        echo '{}' > pipeline_stats.json
        """
}

// Main entry point
workflow {
    // Create input channel
    vcf_ch = Channel.fromPath(params.input)
    
    // Run pipeline
    SNPhylo2(vcf_ch)
    
    // Output results
    SNPhylo2.out.tree.view { "Tree: $it" }
    SNPhylo2.out.report.view { "Report: $it" }
}

"""
Checkpoint and resume functionality for SNPhylo2.

Enables pipeline resumption from intermediate steps.
"""

import json
import hashlib
from pathlib import Path
from typing import Dict, Any, Optional, List
from dataclasses import dataclass, asdict
from datetime import datetime

from snphylo2.utils.logging_utils import get_logger

logger = get_logger()


@dataclass
class StepState:
    """State of a pipeline step."""
    name: str
    status: str  # 'pending', 'running', 'complete', 'failed'
    input_hash: str
    output_files: List[str]
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "StepState":
        return cls(**data)


class CheckpointManager:
    """
    Manages pipeline checkpoints for resumable execution.
    """
    
    def __init__(self, checkpoint_dir: Path, pipeline_id: Optional[str] = None):
        """
        Initialize checkpoint manager.
        
        Args:
            checkpoint_dir: Directory to store checkpoint files
            pipeline_id: Optional unique pipeline identifier
        """
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate pipeline ID from timestamp if not provided
        self.pipeline_id = pipeline_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        self.checkpoint_file = self.checkpoint_dir / f"{self.pipeline_id}.json"
        
        self.state: Dict[str, StepState] = {}
        self._load_checkpoint()
    
    def _load_checkpoint(self) -> None:
        """Load existing checkpoint if available."""
        if self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file, 'r') as f:
                    data = json.load(f)
                
                self.state = {
                    name: StepState.from_dict(state)
                    for name, state in data.get('steps', {}).items()
                }
                
                logger.info(f"Loaded checkpoint: {self.checkpoint_file}")
                logger.info(f"  Completed steps: {self.get_completed_steps()}")
                
            except Exception as e:
                logger.warning(f"Failed to load checkpoint: {e}")
                self.state = {}
    
    def _save_checkpoint(self) -> None:
        """Save current state to checkpoint file."""
        data = {
            'pipeline_id': self.pipeline_id,
            'timestamp': datetime.now().isoformat(),
            'steps': {
                name: state.to_dict()
                for name, state in self.state.items()
            },
        }
        
        # Write to temporary file first
        temp_file = self.checkpoint_file.with_suffix('.tmp')
        with open(temp_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        # Atomic rename
        temp_file.replace(self.checkpoint_file)
    
    def compute_input_hash(self, input_files: List[Path]) -> str:
        """
        Compute hash of input files for change detection.
        
        Args:
            input_files: List of input file paths
            
        Returns:
            MD5 hash of file metadata
        """
        hasher = hashlib.md5()
        
        for file_path in sorted(input_files):
            if file_path.exists():
                # Hash file size and modification time
                stat = file_path.stat()
                hasher.update(f"{file_path}:{stat.st_size}:{stat.st_mtime}".encode())
        
        return hasher.hexdigest()
    
    def is_step_complete(
        self,
        step_name: str,
        input_files: List[Path],
    ) -> bool:
        """
        Check if a step is complete and inputs haven't changed.
        
        Args:
            step_name: Name of the step
            input_files: Current input files
            
        Returns:
            True if step can be skipped
        """
        if step_name not in self.state:
            return False
        
        step = self.state[step_name]
        
        # Check if previously completed
        if step.status != 'complete':
            return False
        
        # Check if inputs changed
        current_hash = self.compute_input_hash(input_files)
        if step.input_hash != current_hash:
            logger.info(f"Step '{step_name}': inputs changed, recomputing")
            return False
        
        # Check if output files exist
        for output_file in step.output_files:
            if not Path(output_file).exists():
                logger.info(f"Step '{step_name}': output missing, recomputing")
                return False
        
        logger.info(f"Step '{step_name}': up to date, skipping")
        return True
    
    def start_step(self, step_name: str, input_files: List[Path]) -> None:
        """
        Mark step as started.
        
        Args:
            step_name: Name of the step
            input_files: Input files for the step
        """
        input_hash = self.compute_input_hash(input_files)
        
        self.state[step_name] = StepState(
            name=step_name,
            status='running',
            input_hash=input_hash,
            output_files=[],
            start_time=datetime.now().isoformat(),
        )
        
        self._save_checkpoint()
        logger.info(f"Started step: {step_name}")
    
    def complete_step(
        self,
        step_name: str,
        output_files: List[Path],
    ) -> None:
        """
        Mark step as complete.
        
        Args:
            step_name: Name of the step
            output_files: Output files produced by the step
        """
        if step_name not in self.state:
            logger.warning(f"Completing unknown step: {step_name}")
            return
        
        self.state[step_name].status = 'complete'
        self.state[step_name].output_files = [str(f) for f in output_files]
        self.state[step_name].end_time = datetime.now().isoformat()
        
        self._save_checkpoint()
        logger.info(f"Completed step: {step_name}")
    
    def fail_step(self, step_name: str, error_message: str) -> None:
        """
        Mark step as failed.
        
        Args:
            step_name: Name of the step
            error_message: Error message
        """
        if step_name not in self.state:
            return
        
        self.state[step_name].status = 'failed'
        self.state[step_name].error_message = error_message
        self.state[step_name].end_time = datetime.now().isoformat()
        
        self._save_checkpoint()
        logger.error(f"Failed step: {step_name} - {error_message}")
    
    def get_completed_steps(self) -> List[str]:
        """Get list of completed step names."""
        return [
            name for name, state in self.state.items()
            if state.status == 'complete'
        ]
    
    def get_failed_steps(self) -> List[str]:
        """Get list of failed step names."""
        return [
            name for name, state in self.state.items()
            if state.status == 'failed'
        ]
    
    def reset_step(self, step_name: str) -> None:
        """
        Reset a step to pending state.
        
        Args:
            step_name: Name of the step to reset
        """
        if step_name in self.state:
            del self.state[step_name]
            self._save_checkpoint()
            logger.info(f"Reset step: {step_name}")
    
    def reset_all(self) -> None:
        """Reset all steps."""
        self.state = {}
        self._save_checkpoint()
        logger.info("Reset all steps")
    
    def get_resume_point(self) -> Optional[str]:
        """
        Determine where to resume pipeline.
        
        Returns:
            Name of first pending step, or None if all complete
        """
        step_order = ['qc', 'filter', 'prune', 'tree', 'report']
        
        for step_name in step_order:
            if step_name not in self.state:
                return step_name
            if self.state[step_name].status != 'complete':
                return step_name
        
        return None


class ResumablePipeline:
    """
    Mixin class for adding resumability to pipelines.
    """
    
    def __init__(self, *args, resume: bool = False, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Initialize checkpoint manager
        checkpoint_dir = Path(self.config.output.directory) / ".checkpoints"
        self.checkpoint = CheckpointManager(
            checkpoint_dir,
            pipeline_id=self.config.output.prefix,
        )
        
        self.resume = resume
        
        if resume:
            resume_point = self.checkpoint.get_resume_point()
            if resume_point:
                logger.info(f"Resuming from step: {resume_point}")
            else:
                logger.info("All steps complete, nothing to resume")
    
    def run_step(self, step_name: str, input_files: List[Path], 
                 output_files: List[Path], func, *args, **kwargs):
        """
        Run a pipeline step with checkpointing.
        
        Args:
            step_name: Name of the step
            input_files: Input files
            output_files: Expected output files
            func: Function to execute
            *args, **kwargs: Arguments for func
            
        Returns:
            Result of func
        """
        # Check if step can be skipped
        if self.resume and self.checkpoint.is_step_complete(step_name, input_files):
            logger.info(f"Skipping completed step: {step_name}")
            return None
        
        # Mark step as started
        self.checkpoint.start_step(step_name, input_files)
        
        try:
            # Execute step
            result = func(*args, **kwargs)
            
            # Verify outputs exist
            missing_outputs = [f for f in output_files if not f.exists()]
            if missing_outputs:
                raise RuntimeError(f"Missing outputs: {missing_outputs}")
            
            # Mark as complete
            self.checkpoint.complete_step(step_name, output_files)
            
            return result
            
        except Exception as e:
            self.checkpoint.fail_step(step_name, str(e))
            raise

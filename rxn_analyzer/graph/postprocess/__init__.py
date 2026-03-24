"""Graph postprocessing helpers.

The new postprocess layer is task-oriented rather than op-oriented.
Current supported mode:
- focus: keep the most important species/reaction nodes and simplify the graph
- context: extract a local neighborhood around one or more seed nodes
- story: extract source-to-target path subgraphs
"""

from .context import run_context_mode
from .focus import run_focus_mode
from .runner import load_postprocess_config, run_postprocess_from_config
from .story import run_story_mode

__all__ = [
    "load_postprocess_config",
    "run_postprocess_from_config",
    "run_focus_mode",
    "run_context_mode",
    "run_story_mode",
]

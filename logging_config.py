import logging
from pathlib import Path


def setup_logging(level=logging.DEBUG):
    """
    Configure logging for the entire project.
    Should be called ONCE from the entry point.
    """

    project_root = Path(__file__).resolve().parent
    log_dir = project_root / "logs"
    log_file = log_dir / "run.log"

    log_dir.mkdir(exist_ok=True)

    # Avoid duplicated handlers (important for Jupyter)
    root = logging.getLogger()
    if root.handlers:
        return

    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    logging.getLogger(__name__).info("Logging initialized")

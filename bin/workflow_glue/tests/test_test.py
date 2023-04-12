"""A dummy test."""

import argparse

from workflow_glue import report_snp


def test():
    """Just showing that we can import using the workflow-glue."""
    assert isinstance(report_snp.argparser(), argparse.ArgumentParser)

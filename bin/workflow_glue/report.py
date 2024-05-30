#!/usr/bin/env python
"""General somatic variation report."""
import json
import os

from dominate.tags import a, code, p
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets.cards import Cards, ICard
from ezcharts.layout.snippets.offcanvas import IOffCanvasClasses

from .util import get_named_logger, wf_parser  # noqa: ABS101


# Reporting function
def populate_report(report, args, **kwargs):
    """Populate the report with the different sections."""
    # Get logger.
    logger = kwargs["logger"]

    # Params.
    params = kwargs["params"]

    # Reports metadata, used to populate the individual cards.
    metadata = {
        "qc": {
            "fn": f"{args.sample_id}.wf-somatic-variation-readQC-report.html",
            "ran": True,
            "exists": None,
            "status": "warning",
            "description": "read alignment statistics",
            "brief": "QC",
        },
        "snv": {
            "fn": f"{args.sample_id}.wf-somatic-snv-report.html",
            "ran": params["snv"],
            "exists": None,
            "status": "fail",
            "description": "small variants calling",
            "brief": "SNV",
        },
        "sv": {
            "fn": f"{args.sample_id}.wf-somatic-sv-report.html",
            "ran": params["sv"],
            "exists": None,
            "status": "fail",
            "description": "structural variant calling",
            "brief": "SV",
        },
        "mod": {
            "fn": f"{args.sample_id}.wf-somatic-mod-report.html",
            "ran": params["mod"],
            "exists": None,
            "status": "fail",
            "description": "modified base analysis",
            "brief": "Mod",
        },
    }

    # Check if the file exists.
    for mode in ["qc", "snv", "sv", "mod"]:
        metadata[mode]["exists"] = os.path.exists(metadata[mode]["fn"])
        if metadata[mode]["exists"] and metadata[mode]["ran"]:
            metadata[mode]["status"] = "pass"
        if not metadata[mode]["exists"] and metadata[mode]["ran"]:
            metadata[mode]["status"] = "warning"
        if not metadata[mode]["exists"] and not metadata[mode]["ran"]:
            metadata[mode]["status"] = "fail"

    # Create list of cards linking to the reports.
    cards = []

    # Create card with link to alignment stats report.
    offcanvas_classes = IOffCanvasClasses()
    offcanvas_classes.offcanvas_button = "btn btn-detail float-end"

    # Create individual cards.
    for mode in ["qc", "snv", "sv", "mod"]:
        # Create link for the data
        option = code(f"--{mode}") if mode != "qc" else "QC"
        link = p("Workflow run without ", option, ".")
        file_link = None
        if metadata[mode]["exists"] and metadata[mode]["ran"]:
            link = a(
                "Report available.",
                href=metadata[mode]["fn"],
            )
        elif not metadata[mode]["exists"] and metadata[mode]["ran"]:
            link = p(
                "Workflow run with ",
                option,
                ", but no report has been generated. ",
                "Check your log for errors.",
            )
        # Create card.
        file_link = {
            "label": "Link",
            "title": "Report for " + metadata[mode]["description"],
            "body": link,
            "classes": offcanvas_classes,
        }
        cards += [
            ICard(
                alias=metadata[mode]["brief"] + " ",
                barcode="",
                body=p(
                    "Report for ", metadata[mode]["description"], " (", option, ")."
                ),
                status=metadata[mode]["status"],
                offcanvas_params=file_link,
            )
        ]

    # Plot the cards.
    logger.info("Start reporting.")
    logger.info("Saving summary section...")
    with report.add_section("Description", "Description"):
        p(
            """
            This report acts as an hub, linking each individual reports generated
            by the workflow.
            Each card links to the reports for each sub-component, if available.
            """
        )
        # Add links to the individual reports in the different cards.
        Cards(columns=len(cards), items=cards)


# Finally, main
def main(args):
    """Run entry point."""
    logger = get_named_logger("report")

    # Instantiate the report.
    report = LabsReport(
        f"{args.sample_id} | Workflow report hub",
        "wf-somatic-variation",
        args.params,
        args.versions,
        args.workflow_version,
        head_resources=[*LAB_head_resources],
    )

    # Load parameters.
    params = json.load(open(args.params))

    # Populate report
    populate_report(
        report=report,
        args=args,
        logger=logger,
        params=params,
    )

    # Save report
    report.write(args.report)
    logger.info(f"Written report to '{args.report}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--sample_id", help="Sample identifier", required=False)
    parser.add_argument(
        "--outdir_path",
        required=True,
        help="Fai index for the reference genome",
    )
    parser.add_argument(
        "--params",
        default=None,
        help="CSV file with workflow parameters",
    )
    parser.add_argument(
        "--versions",
        help="CSV file with software versions",
    )
    parser.add_argument(
        "--workflow_version",
        required=True,
        help="Workflow version",
    )

    return parser

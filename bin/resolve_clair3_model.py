#!/usr/bin/env python
"""Map a basecalling model to a Clair3 model.

An unknown basecalling model or basecalling model without an appropriate
Clair3 model will explode with a large stderr message and exit non-zero.
A happy basecalling model will print a Clair3 model to stdout and exit 0.
"""

# Delegating this to a Python script seems overkill but allows us to
# expand to more complex logic trivially in future.
# Plus I don't want to write this in Groovy right now.

import argparse
import csv
import os
import sys
import textwrap

parser = argparse.ArgumentParser()
parser.add_argument("model_tsv")
parser.add_argument("basecaller_cfg")
parser.add_argument("clair_tool")
parser.add_argument("--liquid_tumor", action='store_true', default=False)
args = parser.parse_args()


def exit_obvious_error(header, error_msg, advice, width=80):
    """Write an obvious looking error to stderr and quit."""
    line = ("-" * width) + '\n'
    msg = (
        f"The input basecaller configuration '{args.basecaller_cfg}' does not "
        f"have a suitable {args.clair_tool} model "
    )
    sys.stderr.write(line)
    sys.stderr.write(f"[CRITICAL ERROR] {header}\n\n")
    sys.stderr.write(textwrap.fill(msg + error_msg, width))
    sys.stderr.write('\n\n')
    sys.stderr.write(textwrap.fill(advice, width))
    sys.stderr.write('\n')
    sys.stderr.write(line)
    sys.exit(os.EX_DATAERR)


with open(args.model_tsv) as tsv:
    for row in csv.DictReader(tsv, delimiter='\t'):
        if row["basecall_model_name"] == args.basecaller_cfg:
            model = row[f"{args.clair_tool.lower()}_model_name"]
            if args.liquid_tumor and row["liquid_model_name_override"] != "-":
                model = row["liquid_model_name_override"]
            reason = row[f"{args.clair_tool.lower()}_nomodel_reason"]
            if model == "-" or model == "":
                # Basecalling model valid but no Clair3 model
                exit_obvious_error(
                    header=f"No appropriate {args.clair_tool} model.",
                    error_msg=f"because {reason}.",
                    advice=(
                        "It is not possible to run the SNP subworkflow "
                        "with this data.\n"
                    )
                )
                break  # exit before here but this keeps my intention obvious
            else:
                # Good model found
                sys.stdout.write(model)
                break
    else:
        # No model found (loop not broken)
        exit_obvious_error(
            header="Unknown basecaller configuration.",
            error_msg=(
                "because the basecaller configuration has not been recognised."
            ),
            advice=(
                "Check your --override_basecaller_cfg has been provided correctly. "
            ),
        )

# CheckAMG/scripts/common/custom_help_formatter.py

import textwrap
import argparse
import shlex

class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    def _fill_text(self, text, width, indent):
        text = textwrap.dedent(text).strip()
        return "\n".join(
            textwrap.fill(line, width, initial_indent=indent, subsequent_indent=indent)
            if line and not line.lstrip().startswith("* ")
            else f"{indent}{line.strip()}"
            for line in text.splitlines()
        )

    def _split_lines(self, text, width):
        lines = []
        for part in text.splitlines():
            if not part.strip():
                lines.append("")
            else:
                lines.extend(textwrap.wrap(part, width))
        return lines
    

def build_rerunnable_command(parser: argparse.ArgumentParser, args: argparse.Namespace) -> str:
    parts = ["checkamg"]
    if getattr(args, "command", None):
        parts.append(str(args.command))

    for action in parser._actions:
        if not action.option_strings:
            continue  # positional
        dest = action.dest
        if dest in {"help", "command", "func", "_raw_argv", "_cli_argv"}:
            continue

        value = getattr(args, dest, None)

        # BooleanOptionalAction: has both --foo and --no-foo
        if isinstance(value, bool) and any(s.startswith("--no-") for s in action.option_strings):
            pos = next((s for s in action.option_strings if s.startswith("--") and not s.startswith("--no-")), None)
            neg = next((s for s in action.option_strings if s.startswith("--no-")), None)
            if pos and neg:
                parts.append(pos if value else neg)
                continue

        # store_true/store_false (single-sided flags)
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            opt = next((s for s in action.option_strings if s.startswith("--")), action.option_strings[0])
            if isinstance(action, argparse._StoreTrueAction):
                if value is True:
                    parts.append(opt)
            else:
                if value is False:
                    parts.append(opt)
            continue

        # lists
        if isinstance(value, (list, tuple)):
            if not value:
                continue
            opt = next((s for s in action.option_strings if s.startswith("--")), action.option_strings[0])
            for v in value:
                if v is None or v == "":
                    continue
                parts.append(opt)
                parts.append(shlex.quote(str(v)))
            continue

        # scalars
        if value is None or value == "":
            continue
        opt = next((s for s in action.option_strings if s.startswith("--")), action.option_strings[0])
        parts.append(opt)
        parts.append(shlex.quote(str(value)))

    return " ".join(parts)

from pathlib import Path
import re
import sys

src = Path(sys.argv[1])
dst = Path(sys.argv[2])

name = sys.argv[3]
tend = sys.argv[4]
restart = sys.argv[5]
write_restart = sys.argv[6].lower() == "true"
restart_out = sys.argv[7]

text = src.read_text()

def set_ini(text, section, key, value):

    lines = text.splitlines()

    sec_re = re.compile(
        rf"^\s*\[{re.escape(section)}\]\s*$",
        re.I
    )

    key_re = re.compile(
        rf"^\s*{re.escape(key)}\s*=",
        re.I
    )

    sec_idx = None

    for i, line in enumerate(lines):
        if sec_re.match(line):
            sec_idx = i
            break

    if sec_idx is None:
        if lines and lines[-1].strip():
            lines.append("")
        lines.append(f"[{section}]")
        lines.append(f"{key} = {value}")
        return "\n".join(lines) + "\n"

    end = len(lines)

    for j in range(sec_idx + 1, len(lines)):
        if re.match(r"^\s*\[.+\]\s*$", lines[j]):
            end = j
            break

    for j in range(sec_idx + 1, end):
        if key_re.match(lines[j]):
            lines[j] = f"{key} = {value}"
            return "\n".join(lines) + "\n"

    lines.insert(end, f"{key} = {value}")

    return "\n".join(lines) + "\n"


text = set_ini(
    text,
    "TimeLoop",
    "TEnd",
    tend
)

text = set_ini(
    text,
    "Problem",
    "Name",
    name
)

text = set_ini(
    text,
    "Restart",
    "Enable",
    "true"
)

text = set_ini(
    text,
    "Restart",
    "File",
    restart
)

text = set_ini(
    text,
    "Restart",
    "WriteAtEnd",
    "true" if write_restart else "false"
)

if write_restart:
    text = set_ini(
        text,
        "Restart",
        "OutputFile",
        restart_out
    )

dst.write_text(text)

print("Wrote:", dst)

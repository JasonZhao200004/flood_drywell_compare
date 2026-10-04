#!/usr/bin/env bash

OUT="$1"
RUN="$2"
PY="$3"

PIDFILE="$OUT/${RUN}.pid"

F100="$OUT/${RUN}-00100.vtu"
F103="$OUT/${RUN}-00103.vtu"

echo "103-h atmospheric-basin watchdog started: $(date)"

while [[ ! -f "$F103" ]]; do

    if [[ -f "$PIDFILE" ]]; then
        PID="$(cat "$PIDFILE")"

        if ! kill -0 "$PID" 2>/dev/null; then
            echo "MODEL STOPPED BEFORE 103 h"
            exit 1
        fi
    fi

    sleep 60
done

"$PY" - "$F100" "$F103" <<'PY'
import sys,re
import numpy as np
import meshio


def norm(s):
    return re.sub(
        r"[^a-z0-9]",
        "",
        s.lower()
    )


def get(path):

    m=meshio.read(path)

    tri=[
        i for i,b in enumerate(m.cells)
        if b.type=="triangle"
    ]

    for name,parts in m.cell_data.items():

        n=norm(name)

        if "o2tag" not in n or "gas" not in n:
            continue

        x=[]

        for i in tri:
            x.append(
                np.asarray(
                    parts[i],
                    float
                ).reshape(-1)
            )

        return np.concatenate(x)

    raise RuntimeError(
        "x^O2tag_gas not found"
    )


a=get(sys.argv[1])
b=get(sys.argv[2])

m100=float(np.max(np.abs(a)))
m103=float(np.max(np.abs(b)))

n103=int(np.sum(b>1e-12))

print("============================================================")
print("ATMOSPHERIC BASIN TAG CHECK")
print("100 h max|gas tag| =",f"{m100:.12e}")
print("103 h max|gas tag| =",f"{m103:.12e}")
print("103 h cells >1e-12 =",n103)

if m100>1e-10:
    raise SystemExit(
        "FAIL: tag exists before drainage atmosphere opens"
    )

if m103<=1e-10 or n103==0:
    raise SystemExit(
        "FAIL: basin atmosphere/tag did not activate"
    )

print("103-h CHECK: PASS")
print("Atmospheric O2 is entering through the exposed basin.")
print("============================================================")
PY

RC=$?

if [[ "$RC" -ne 0 ]]; then

    echo "TAG CHECK FAILED — stopping formal run."

    PID="$(cat "$PIDFILE" 2>/dev/null || true)"

    if [[ -n "$PID" ]] && kill -0 "$PID" 2>/dev/null; then
        kill "$PID"
    fi

    exit "$RC"
fi

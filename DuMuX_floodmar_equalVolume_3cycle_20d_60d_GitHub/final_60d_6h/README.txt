FINAL FLOOD MAR SAME-VOLUME DATASET

Simulation:
    0 to 1440 h = 60 days

Final schedule:
    C1 shutdown start = 96.000000000000 h
    C2 shutdown start = 576.164044805955 h
    C3 shutdown start = 1060.953978421305 h

Target injection:
    3956.009528 m3/cycle

Archive:
    nominal interval = 6 h
    states = 241

Sources:
    0-576 h:
        existing verified final-compatible trajectory

    582-954 h:
        Mac continuation from restart 512.061111111 h
        using final C2 schedule

    960-1026 h:
        existing trajectory generated after final C2
        and before final C3 shutdown

    1032-1440 h:
        Mac continuation from restart 1030 h
        using final C3 schedule

IMPORTANT:
    PVD timesteps preserve ACTUAL simulation times.
    The 512.061111-h restart causes the Part-1 states
    to differ from nominal integer-hour targets by about
    3.667 minutes.

    See OUTPUT_TIME_MAPPING.csv for exact times.

Maximum nominal/actual offset:
    8.500000 min

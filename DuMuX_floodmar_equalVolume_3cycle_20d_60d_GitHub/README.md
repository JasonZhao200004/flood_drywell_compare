# Flood MAR final same-volume continuation bundle

Final target = 3956.009528 m3/cycle

Final shutdown starts:
- C1 = 96.000000000000 h
- C2 = 576.164044805955 h
- C3 = 1060.953978421305 h

Already verified and packaged:
- 0–576 h, every 6 h: 97 states
- 960–1026 h, every 6 h: 12 states

Mac must generate:
- 582–954 h from the 512.061111 h restart
- 1032–1440 h from the 1030 h restart

After merging:
- 0–1440 h
- every 6 h
- 241 states

The old C1C2formal trajectory is used ONLY through 576 h, before
the final C2 schedule diverges from the old schedule.

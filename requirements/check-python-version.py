#!/usr/bin/env python
import sys
if not (sys.version_info[0] == 3 and sys.version_info[1] == 6):
    print(
        '\033[31mERROR\033[0m: '
        'Pinned dependencies should only be set from Python 3.6',
        file=sys.stderr)
    sys.exit(1)

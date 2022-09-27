#!/usr/bin/env python3
#
# Builds the documentation (in a subprocess). Can be used as a test.
#
import inspect
import os
import subprocess
import sys


def doctest_sphinx(root=os.curdir):
    """Runs sphinx-build in a subprocess, checking that it can be invoked without producing errors."""

    print('Checking if docs can be built.')
    sys.stdout.flush()  # flush output to see the printed output _before_ the subprocess output is shown

    # Build the HMTL docs in a subprocess
    p = subprocess.Popen([
        'sphinx-build',
        '-b',
        'html',
        os.path.join(root, 'source'),
        os.path.join(root, 'build', 'html'),
        '-W',
    ])

    # Wait for process to terminate, and check status
    try:
        ret = p.wait()
    except KeyboardInterrupt:
        # Allow Ctrl-C to interrupt
        try:
            p.terminate()
        except OSError:
            pass
        p.wait()
        print('Build CANCELLED')
        sys.exit(1)

    # Check returned status
    if ret != 0:
        print('Build FAILED')
        sys.exit(ret)
    print('Build OK')


if __name__ == '__main__':
    # Get path to docs directory (so that this can be called from the root directory).
    try:
        frame = inspect.currentframe()
        path = os.path.abspath(os.path.dirname(inspect.getfile(frame)))
    finally:
        # Always manually delete frame: https://docs.python.org/2/library/inspect.html#the-interpreter-stack
        del frame
    path = os.path.relpath(path, os.getcwd())

    # Run
    doctest_sphinx(path)

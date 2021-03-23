#!/usr/bin/env python3
#
# Builds the documentation (in a subprocess). Can be used as a test.
#
import inspect
import os
import subprocess
import sys


def doctest_sphinx(root=os.curdir):
    """
    Runs sphinx-build in a subprocess, checking that it can be invoked without
    producing errors.
    """
    print('Checking if docs can be built.')
    sys.stdout.flush()
    p = subprocess.Popen([
        'sphinx-build',
        '-b',
        'html',
        os.path.join(root, 'source'),
        os.path.join(root, 'build', 'html'),
        '-W',
    ])
    try:
        ret = p.wait()
    except KeyboardInterrupt:
        try:
            p.terminate()
        except OSError:
            pass
        p.wait()
        print('Build OK')
        sys.exit(1)
    if ret != 0:
        print('Build FAILED')
        sys.exit(ret)


if __name__ == '__main__':
    # Get path to docs directory (so that this can be called from the root
    # directory).
    try:
        frame = inspect.currentframe()
        path = os.path.abspath(os.path.dirname(inspect.getfile(frame)))
    finally:
        # Always manually delete frame
        # https://docs.python.org/2/library/inspect.html#the-interpreter-stack
        del(frame)
    path = os.path.relpath(path, os.getcwd())

    # Run
    doctest_sphinx(path)

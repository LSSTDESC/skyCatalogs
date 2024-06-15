"""
Utilities for top-level scripts
"""
from datetime import datetime as dt
import sys
import logging

__all__ = ['print_callinfo', 'log_callinfo', 'callinfo_to_dict', 'print_date',
           'print_dated_msg', 'TIME_TO_SECOND_FMT']

TIME_TO_SECOND_FMT = '%Y-%m-%d %H:%M:%S'

def print_callinfo(prog, args):
    """
    Print information about how a script using argparse was called

    Parameters
    ----------
    prog   program name, typically sys.argv[0]
    args   object returned by ArgumentParser.parse_args()
    """

    print('{}   {}  invoked  with arguments'.format(dt.now().strftime(TIME_TO_SECOND_FMT), prog))
    for e in dir(args):
        if not e.startswith('_'):
            nm = 'args.' + e
            print('{}: {}'.format(e, eval(nm)))

    sys.stdout.flush()

def log_callinfo(prog, args, logname):
    """
    Write information about how a script using argparse was called to log

    Parameters
    ----------
    prog   program name, typically sys.argv[0]
    args   object returned by ArgumentParser.parse_args()
    logname string
    """

    logger = logging.getLogger(logname)
    log_out = '{}  invoked  with arguments\n'.format(prog)
    for e in dir(args):
        if not e.startswith('_'):
            nm = 'args.' + e
            log_out += '    {}: {}\n'.format(e, eval(nm))
    logger.info(log_out)

def callinfo_to_dict(args):
    """
    Make a dict out of program arguments.  Each option value is
    either a simple atomic type or a list
    """
    a_dict = dict()
    for e in dir(args):
        if not e.startswith('_'):
            v = eval('args.' + e)
            a_dict[e] = v
    return a_dict

def print_date(to_second=True, file=None):
    """
    Print current time (by default only to nearest second) and flush output
    Print to file if supplied, else to std.out
    """
    if to_second:
        print(dt.now().strftime(TIME_TO_SECOND_FMT), file=file, flush=True)
    else:
        print(dt.now(), file=file, flush=True)
def print_dated_msg(msg, to_second=True, file=None):
    """
    Print current time (by default only to nearest second) and flush output
    Print to file if supplied, else to std.out
    """
    if to_second:
        print(dt.now().strftime(TIME_TO_SECOND_FMT), ' ', msg,  file=file, flush=True)
    else:
        print(dt.now(), ' ', msg, file=file, flush=True)

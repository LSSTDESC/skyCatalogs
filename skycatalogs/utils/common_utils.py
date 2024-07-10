"""
Utilities for top-level scripts
"""
from datetime import datetime as dt
import sys
import logging

__all__ = ['log_callinfo', 'callinfo_to_dict', 'print_date',
           'print_dated_msg', 'TIME_TO_SECOND_FMT']

TIME_TO_SECOND_FMT = '%Y-%m-%d %H:%M:%S'


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
    for k,v in dict(sorted(args._get_kwargs())).items():
        log_out += '    {}: {}\n'.format(k, v)
    logger.info(log_out)

def callinfo_to_dict(args):
    """
    Make a dict out of program arguments.  Each option value is
    either a simple atomic type or a list
    """
    return dict(args._get_kwargs())

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

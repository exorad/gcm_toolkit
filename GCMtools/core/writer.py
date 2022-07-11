# ==============================================================
#                       Writer Subroutines
# ==============================================================
#  This file contains all functions to handle informative and
#  diagnostic output during the run. If writer.on = False, no
#  output will be produced and GCMtools will run silently. If a
#  file_name is given, it will save the output to this file.
# ==============================================================

class _writer:
    """
    This class saves all properties of the writing routines. In
    particular it remembers if output should be given to files
    or the console, the colors of the print and the layout.
    """

    # Set whether output is muted or not
    on = False

    # set wheter output is to console (None) or to file
    file_name = None

    # color settings
    ENDC = '\033[0m'
    DEFAULT = '\033[94m'
    ERROR = '\033[91m'
    WARN = '\033[93m'

    # spare colors and layouts if needed
    # HEADER = '\033[95m'
    # OKCYAN = '\033[96m'
    # OKGREEN = '\033[92m'
    # BOLD = '\033[1m'
    # UNDERLINE = '\033[4m'

    # max line length of print
    line_length = 80

    # max length before words get cut to enforce line length
    line_maxoff = 10

    # spacing for 1 unite of indent
    spacer = "   "


def writer_setup(typ):
    """
    Set-Up class for all writer function. Sets the correct
    properties in the _writer class. Careful: if called mid
    run, it overwrites the log file.

    Parameters
    ----------
    typ : str
        Type of output, can either be:
        'off': no output
        'on': output to console
        str: file path to where log should be saved
    """

    # mute all outputs
    if typ == 'off':
        _writer.on = False
        _writer.file_name = None

    # enable all outputs to console
    elif typ == 'on':
        _writer.on = True
        _writer.file_name = None

    # enable all outputs to file
    else:
        _writer.on = True
        _writer.file_name = typ
        f = open(typ, 'w')
        f.close()


def write_status(tag, message):
    """
    Write a status output that is auto layouted according to
    _writer set up. Careful: an ERROR status will call a
    ValueError.

    Parameters
    ----------
    tag : str
        tag of the status. Some have pre-defined meaning:
        'STAT', 'INFO', 'E-INFO', 'WARN', 'ERROR'
    message : str
        The message to be printed.
    """

    # check if verbosity is whished
    if not _writer.on:
        return

    # define colors and layout according to tags
    color = _writer.DEFAULT
    space = ""
    if tag == 'STAT':
        color = _writer.DEFAULT
    if tag == 'INFO':
        color = _writer.DEFAULT
        space = _writer.spacer
    if tag == 'E-INFO':
        color = _writer.ERROR
    if tag == 'WARN':
        color = _writer.WARN
    if tag == 'ERROR':
        color = _writer.ERROR

    # add line breaks after _writer.line_length characters
    line = space + "[" + tag + "] " + message
    tag_length = ' ' * (len(space + "[" + tag + "] "))
    liner = ""
    no_loop = True
    while len(line) > _writer.line_length:
        no_loop = False
        for i in range(_writer.line_maxoff):
            if line[_writer.line_length - i] == " ":
                liner += line[:_writer.line_length - i + 1] + "\n"
                line = tag_length + line[_writer.line_length - i + 1:]
                break
        else:
            liner += line[:_writer.line_length + 1] + "\n"
            line = tag_length + line[_writer.line_length + 1:]
    liner += line

    # write message to terminal
    if _writer.file_name == None:
        print(color + liner + _writer.ENDC)

    # write message to file
    else:
        f = open(_writer.file_name, 'a')
        f.write(liner + "\n")
        f.close()

    if tag == 'ERROR':
        raise ValueError(liner)


def write_message(message, color=None, spacing=0):
    """
    Write a status output that is auto layouted according to
    _writer set up. Careful: an ERROR status will call a
    ValueError.

    Parameters
    ----------
    message : str
        The message to be printed.
    color : str
        Color of the message, only options availble are:
        'WARN', 'ERROR'
        Otherwise, the default color will be used
    spacing : int
        left spacing of the message in number of whtiespaces
    """

    # check if verbosity is whished
    if not _writer.on:
        return

    # set color scheme
    col = _writer.DEFAULT
    if color == 'WARN':
        col = _writer.WARN
    elif color == 'ERROR':
        col = _writer.ERROR

    # add line breaks after _writer.line_length characters
    tag_length = ' ' * spacing
    line = tag_length + message
    liner = ""
    no_loop = True
    while len(line) > _writer.line_length:
        no_loop = False
        for i in range(_writer.line_maxoff):
            if line[_writer.line_length - i] == " ":
                liner += line[:_writer.line_length - i + 1] + "\n"
                line = tag_length + line[_writer.line_length - i + 1:]
                break
        else:
            liner += line[:_writer.line_length + 1] + "\n"
            line = tag_length + line[_writer.line_length + 1:]
    liner += line

    # write message to terminal
    if _writer.file_name == None:
        print(col + liner + _writer.ENDC)

    # write message to file
    else:
        f = open(_writer.file_name, 'a')
        f.write(liner + "\n")
        f.close()


def write_hline(character='=', color=None):
    """
    Write a single horizontal line.

    Parameters
    ----------
    character : str
        line character (If given a string, only the first
        character will be used.)
    color : str
        Color of the message, only options available are:
        'WARN', 'ERROR'
        Otherwise, the default color will be used.
    """

    # check if verbosity is whished
    if not _writer.on:
        return

    # set color scheme
    col = _writer.DEFAULT
    if color == 'WARN':
        col = _writer.WARN
    elif color == 'ERROR':
        col = _writer.ERROR

    # define output
    liner = character[0] * _writer.line_length

    # write message to terminal
    if _writer.file_name == None:
        print(col + liner + _writer.ENDC)

    # write message to file
    else:
        f = open(_writer.file_name, 'a')
        f.write(liner + "\n")
        f.close()

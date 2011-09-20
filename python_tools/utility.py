
def line_break():
    """Print a line break"""
    
    print "===================="

def mkdir(dir_name):
    """Make a directory but if it exists then ignore the error.
    c.f. http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    """
    try:
        os.__name__
        errno.__name__
    except NameError:
        import os, errno
    
    try:
        os.makedirs('a/b/c')
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise


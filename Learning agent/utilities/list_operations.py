from itertools import zip_longest


def grouper(iterable, n, fillvalue=None):
    """Collects data into fixed-length chunks or blocks.

    Taken from "https://docs.python.org/3/library/itertools.html.

    Parameters
    ------------
    iterable:     the object that will be iterated through in n-sized chunks.
    n:            (int) size of the chunk.
    fillvalue:    (optional) the value appended to the iterable if the iterable
                  is not divisible by n.

    Returns
    ---------
    (iterable) an iterable where each item is of size n.
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

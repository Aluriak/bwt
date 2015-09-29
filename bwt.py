"""
Burrows-Wheeler Transform pure Python implementation.

Here the bwt keyword references the Burrows-Wheeler Transform,
obtained from a text (generally named 'text').
Many metadata, like ranks, alphabet and counts of char, are generated
and used by different algorithms implemented in this module.
Save these metadata and give them to the API functions can lead to important
efficacity gains.

Main API is the following functions:
    - bwt, that gives the bwt string of the given text.
    - construct_text, that generate the text from the bwt.
    - find, that give the number of occurences of a query in the text,
        or yield the postions of these occurences if enough metadata provided.

"""
from operator import itemgetter
from decorator import decorator
from functools import partial
from collections import Counter


LAST_CHARACTER = '\0'  # used as tail-end character in BWT algorithms


def validate_tail(char=LAST_CHARACTER):
    """Decorator that add assure the presence of given char
    at the end of the text keyword argument

    """
    @decorator
    def _validate_tail(func, *args, **kwargs):
        """Decorator effectively returned as decorator of func

        """
        # note that kwargs will always be empty of positionnal arguments,
        # thanks to the decorator.decorator decorator.
        # this same decorator automatically wraps the decoratorated function.
        args = list(args)
        if args[0][-1] != char:
            args[0] += char
        return func(*args, **kwargs)
    return _validate_tail


@validate_tail()
def generate_suffix_array(text):
    """Return the suffix array of text

    >>> suffix_array('TAATATA')
    [(7, '\x00'), (6, 'A\x00'), (1, 'AATATA\x00'), (4, 'ATA\x00'), (2, 'ATATA\x00'), (5, 'TA\x00'), (0, 'TAATATA\x00'), (3, 'TATA\x00')]

    """
    return tuple(idx for idx, suffix in sorted(
        ((i, text[i:]) for i in range(len(text))),
        key=itemgetter(1)
    ))


@validate_tail()
def transform(text, suffix_array=None):
    """Yields the Burrows-Wheeler Transformed text

    if suffix_array is not given (set to None),
    it will be generated before any treatment.

    >>> tuple(bwt.generate('ACAACG'))
    ('G', 'C', '\x00', 'A', 'A', 'A', 'C')
    >>> tuple(bwt.generate('CATGTATT'))
    ('T', 'C', 'T', '\x00', 'T', 'T', 'G', 'A', 'A')

    """
    if suffix_array is None:
        suffixes_order = generate_suffix_array(text)
    return (
        text[suffixes_order[it] - 1] if suffixes_order[it] > 0 else '\0'
        for it in range(len(text))
    )


def generate_ranks(bwt):
    """Yields ranks of chars in given bwt, and finally the counter

    >>> tuple(bwt.generate_ranks(
        'TCT\0TTGAA')
    )
    (1, 1, 2, 1, 3, 4, 1, 1, 2)
    >>> tuple(bwt.generate_ranks(bwt.bwt('CATGTATT')))
    (1, 1, 2, 1, 3, 4, 1, 1, 2)

    """
    count = Counter()
    for char in bwt:
        count[char] += 1
        yield count[char]
    yield count


def left_forward(char, k, counts, alphabet=None):
    """Return the index of the k-th char in the list of all char
    of the studied text, ordered as in given alphabet or lexicographical order.

    Alphabet can be interpolated from counts,
    that is a dict character:count in text, by use the lexicographical order.

    >>> bwt.left_forward('C', 3, counts={'A':100, 'B':10, 'C':5, 'D':1000})
    112
    # in this example, the alphabet is 'ABCD'
    # the third C is after the A, the B, and the k-1 firsts C : 100 + 10 + 3-1

    """
    if alphabet is None:
        alphabet = sorted(counts.keys())
    # get all letters that are lexicographically smaller than char
    previous_letters = alphabet[:alphabet.index(char)]
    # sum count of these letters in the text
    offset = sum(counts[letter] for letter in previous_letters)
    # the k-th character is in the list of characters,
    # k-1 places after the last smaller letter
    return offset + k - 1


def generate_metadata(bwt, ranks=None, counts=None, alphabet=None):
    """Return ranks, counts and alphabet given or deduced from given bwt"""
    if alphabet is None:
        alphabet = tuple(sorted(tuple(set(bwt))))
    if ranks is None or counts is None:
        *ranks, counts = tuple(generate_ranks(bwt))
    # print('ALPHABET, RANKS, COUNTS:', alphabet, ranks, counts)
    return ranks, counts, alphabet


def construct_text(bwt, *, ranks=None, counts=None, alphabet=None):
    """Yields the text described by the bwt.

    ranks, counts and alphabet can be given. If not, they are generated
    before any other treatment.

    """
    ranks, counts, alphabet = generate_metadata(bwt, ranks, counts, alphabet)
    next_line = partial(left_forward, alphabet=alphabet, counts=counts)
    line = 0
    text = bwt[line]
    while bwt[line] != '\0':
        line = next_line(bwt[line], ranks[line])
        text = bwt[line] + text
    if text[0] == LAST_CHARACTER:
        # avoid using the LAST_CHARACTER (marking the tail end)
        # (NB: reconstruction place it at the beginning)
        return text[1:]
    else:
        return text


def find(query, bwt, suffix_array=None, *, ranks=None, counts=None, alphabet=None):
    """Return number of occurences of query in the text described by given bwt

    Use the Ferragina-Mozzina Index method.
    If suffix_array, the suffix array of the initial text, is given,
    this function turns into a generator that yields positions
    of the query in the text.

    """
    # data definition
    ranks, counts, alphabet = generate_metadata(bwt, ranks, counts, alphabet)
    next_line = partial(left_forward, alphabet=alphabet, counts=counts)
    # init
    start, stop = 0, len(bwt) - 1

    for target in reversed(query):
        # find the minimal possible char in bwt
        idmin = start
        while bwt[idmin] != target and idmin <= stop: idmin += 1
        # same for the maximal one
        idmax = stop
        while bwt[idmax] != target and idmax >= start: idmax -= 1
        # is something found ?
        if idmax - idmin >= 0:
            start = next_line(target, ranks[idmin])
            stop  = next_line(target, ranks[idmax])
            if suffix_array is not None:
                pass  # TODO: implement the positions yielding
        else:  # nothing found: stop the search here
            break
    # return the number of occurences only if all occurences are not yet yielded
    if suffix_array is None:
        return stop - start + 1


if __name__ == '__main__':

    TXT = 'CATGTATT'
    BWT = tuple(transform(TXT))
    assert TXT == construct_text(BWT)
    print('Text:', construct_text(BWT))
    print('BWT :', BWT)
    for query in ('T', 'AT', 'CAT'):
        print(query.rjust(4) + ':', find(query, BWT), 'occurences')


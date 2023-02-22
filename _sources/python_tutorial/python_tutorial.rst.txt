 .. _python-tutorial-label:

===============
Python Tutorial
===============

The python tutorial of 'start science here!' certainly is the longest and most exhaustive part. Python is for a computer scientist (computational chemists to be exact) the toolkit with which they do their work. Just like a carpenter needs their hammers or an organic chemist needs their flasks, we use python for our work.

First things first: You don't have to work through the python tutorial from start to finish. For a start the following Quick Reference might suffice. But if you find the time you can head over to binder and start with the notebooks in the `python_tutorial` folder. You can stop at any point and do some other tutorials or work on your own project. But let me tell you, that I learned new stuff while composing the tutorial and I was using python daily for 3 years at this point.

Introduction to Jupyter Notebooks
=================================

Before clicking on that binder link and starting a python session read this short introduction about jupyter notebooks, how to work with them and how to execute code inside of them.

.. toctree::

   introduction.nblink

Interactive Tutorial Notebooks
==============================

Available on binder:

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/kevinsawade/start-science-here/HEAD?urlpath=%2Ftree%2Ffilepath=python_tutorial
   :alt: Binder link to repo.

Non-interactive static notebooks
================================

An offline, non-interactive version of the notebooks can be found via this menu:

.. toctree::
   :maxdepth: 1
   :caption: Offline notebooks

   basics_00_datatypes.nblink
   basics_01_functions_classes.nblink
   intermediate_02_OOP.nblink

Quick Reference
===============

The python quick reference or cheat sheet was taken from https://github.com/justmarkham/python-reference and adjusted for start science here. The remainder of this page will contain this quick reference. Each topic has a direct link to the examples in binder and an offline version of the examples.

Imports
-------

Import modules with or without alias

.. code-block:: python

   import math
   import numpy as np
   print(math.pi)
   print(np.pi)

Multiple module imports and multiple imports from a module

.. code-block:: python

   import os, sys
   from numpy.random import random, randint
   print(os.getcwd())
   print(sys.version)
   print(random())
   print(randint(100))

Import all objects from a module (this is generally discouraged, as it makes your namespace messy)

.. code-block:: python

   from math import *
   print(pi)
   print(__doc__) # example of messy namespace

Check the attributes and methods of any object.

.. code-block:: python

   import math
   from numpy import ndarray
   print(dir(math))
   print(dir(ndarray))

Data Types
----------

To quickly start into an interactive environment with infos about data types:

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/kevinsawade/start-science-here/HEAD?urlpath=%2Ftree%2F?filepath=python_tutorial%2Fquick_reference%2Fqr_01_data_types.ipynb
   :alt: Binder link to qr_01_data_types.ipynb

.. image:: https://img.shields.io/badge/offline-notebook-orange
   :target: qr_01_data_types.html
   :alt: Link to offline 01_data_types.ipynb

Determine the type of an object:

.. code-block:: python

   import math
   print(type(math))    # module
   from numpy import ndarray
   print(type(ndarray)) # classes will print `type`
   from numpy.random import random
   print(type(random))  # builtin_function_or_method
   print(type(1))       # int (integer)
   print(type(1.1))     # float (floating point number)
   print(type('hello')) # str (string)
   print(type(True))    # bool (boolean, True/False values)
   print(type(None))    # NoneType (NoneType is a singleton)

Check the type of an object:

.. code-block:: python

   isinstance(2, int) # True
   isinstance(2.2, (float, int)) # True

Converting between types:

.. code-block:: python

   float(2) # 2.0
   int(2.9) # 2
   int(3.9) # 3
   str(2.2) # '2.2'


Boolean values of built-in types:

.. code-block:: python

   bool(0)    # False
   bool(None) # False
   bool('')   # False
   bool([])   # False (empty list)
   bool({})   # False (empty dict)

True for variables with values and non-empty containers:


.. code-block:: python

   bool(2)     # True
   bool([2])   # True
   bool('two') # True

Math
----

Arithmetic operations

.. code-block:: python

   10 + 4        # addition
   10 - 4        # subtraction
   10 * 4        # multiplication
   10 / 4        # division (python 2.7 has a quirk here and wont return the expected result)
   10 ** 4       # exponent
   10 // 4       # floor division
   10 % 4        # modulo
   divmod(10, 4) # returns tuple of floor division and modulo

Comparsions and boolean operations
----------------------------------

Assigning objects to variables:

.. code-block:: python

   a = 5

Arithmetic comparisons:

.. code-block:: python

   x > 3  # True
   x >= 3 # True
   x != 3 # True
   x == 3 # False

Difference between ``==`` and ``is``:

.. code-block:: python

   a = []
   b = []
   a == b # True
   a is b # False
   a is a # True

Boolean operations and combinations:

.. code-block:: python

   5 > 3 and 6 > 3 # True
   5 > 3 or 5 < 3  # True
   not False       # True
   False or not False and True # ?

Control Flow
------------

If/elif/else:

.. code-block:: python

   if x > 0:
       print('positive')
   elif x == 0:
       print('zero')
   else:
       print('negative')

For loops (can contain ``else`` to do something after the loop is complete):

.. code-block:: python

   for i in range(5):
       print(i ** 2)
   else:
       print("for loop closed")

Continue and break:

.. code-block:: python

   for i in range(10):
       if i == 0:
           # don't print 0
           continue
       if i == 5:
           # break loop after 5
           break
       print(i)

While loops continue running until either ``break`` occurs or the statement becomes false.

.. code-block:: python

   max_square = 50
   i = 0
   while i ** 2 <= max_square:
       print(i ** 2)
       i += 1
       if i > 100:
           break # security break

Lists
-----

**Properties of lists: ordered, iterable, mutable**

Empty lists (``bool([]) = False``):

.. code-block:: python

   empty_list = []
   empty_list = list()


Create a list:

.. code-block:: python

   simpsons = ['homer', 'marge', 'bart']

Examine a list:

Access an element in the list with indexing (``[index]``) and multiple elements with slicing (``[start:stop:step]``). Get length with the built-in ``len()`` function.

.. code-block:: python

   print(simpsons[0])
   print(len(simpsons))

Modify lists (does not return lists, but rather changes the list in-place).

.. code-block:: python

   simpsons.append('lisa')                # add a single element
   simpsons.extend(['itchy', 'scratchy']) # add all elements from another list
   simpsons.insert(0, 'maggie')           # add element at index
   simpsons.remove('bart')                # remove this specific element
   simpsons.pop(0)                        # remove element at index and return that element
   del simpsons[0]                        # remove that element and not return it
   simpsons[0] = 'krusty'                 # overwrite indexed element


Addition of lists returns new lists:

.. code-block:: python

   neighbors = simpsons + ['ned', 'rod', 'todd']

Count element in list.

.. code-block:: python

   simpsons.count('lisa')      # counts the number of instances
   simpsons.index('itchy')     # returns index of first instance

List slicing:


.. code-block:: python

   weekdays = ['mon', 'tues', 'wed', 'thurs', 'fri']
   weekdays[0]         # element 0
   weekdays[0:3]       # elements 0, 1, 2
   weekdays[:3]        # elements 0, 1, 2
   weekdays[3:]        # elements 3, 4
   weekdays[-1]        # last element (element 4)
   weekdays[::2]       # every 2nd element (0, 2, 4)
   weekdays[::-1]      # backwards (4, 3, 2, 1, 0)

Alternative to ``[::-1]``: list(reversed(...))

.. code-block:: python

   list(reversed(weekdays))

Sorting lists:

Sorts a list (alphanumerical) in-place. This overwrites the original list and changes the order of the elements:

.. code-block:: python

   simpsons.sort()
   simpsons.sort(reverse=True)     # sort in reverse
   simpsons.sort(key=len)          # sort by a key

To get a *new* sorted list, without changing the original:

.. code-block:: python

   sorted(simpsons)
   sorted(simpsons, reverse=True)
   sorted(simpsons, key=len)

Insert into already sorted list:

.. code-block:: python

   num = [10, 20, 40, 50]
   from bisect import insort
   insort(num, 30)

Creating a variable from a new list **does not** copy it. Instead you are creating a reference to the same list:

.. code-block:: python

   same_num = num
   print(same_num)
   same_mun[0] = 0
   print(num)
   print(same_num)

If you want to keep ``num`` and copy it to ``same_num``, here are three ways:

.. code-block:: python

   same_num = num[:]
   same_num = list(num)
   import copy
   same_num = copy.deepcopy(num)

Copying elements more thorough:

.. code-block:: python

   import copy

   a = [1, 2, 3]
   b = [4, 5, 6]
   c = [a, b]

Normal assignment:

.. code-block:: python

   d = c

   print(c == d)          # True - equality == is not identity is
   print(c is d)          # True - d is the same object as c
   print(c[0] is d[0])    # True - d[0] is the same object as c[0]

Shallow (normal) copy constructs a new compound object but keeps the underlying references to already existing objects:

.. code-block:: python

   d = copy.copy(c)

   print(c == d)          # True - equality == is not identity is
   print(c is d)          # False - d is now a new object
   print(c[0] is d[0])    # True - d[0] is the same object as c[0]

Deepcopy constructs a new compound object and then, recursively, inserts copies into the new object found in the original object:

.. code-block:: python

   d = copy.deepcopy(c)

   print(c == d)          # True - equality == is not identity is
   print(c is d)          # False - d is now a new object
   print(c[0] is d[0])    # False - d[0] is now a new object


Tuples
------

**Properties of lists: ordered, iterable, immutable**

The main difference between tuples and lists is, that tuples are immutable. You can not change an element inside a tuple with assignment ``(1, 2, 3)[2] = 4``.

Empty tuples:

.. code-block:: python

   a = tuple()
   b = ()
   print(bool(a), bool(b))

Tuples with values:

.. code-block:: python

   digits = (0, 1, 'two')          # create a tuple directly
   digits = tuple([0, 1, 'two'])   # create a tuple from a list
   zero = (0,)                     # trailing comma is required to indicate it's a tuple

Getting values from a tuple:

.. code-block::

   digits[2]           # returns 'two'
   len(digits)         # returns 3
   digits.count(0)     # counts the number of instances of that value (1)
   digits.index(1)     # returns the index of the first instance of that value (1)

Assigning values fails:

.. code-block::

   digits[2] = 2       # throws an error

You need to compose new tuples to make that work:

.. code-block:: python

   new_digits = (*digits[:1], 2) # more on the asterisk in Packing and Unpacking

Or concatenate tuples:

.. code-block:: python

   digits = digits + (3, 4)

Multiplication (also works with lists):

.. code-block:: python

   (3, 4) * 2          # returns (3, 4, 3, 4)

Sort a list of tuples:

.. code-block:: python

   tens = [(20, 60), (10, 40), (20, 30)]
   sorted(tens)        # sorts by first element in tuple, then second element
                       #   returns [(10, 40), (20, 30), (20, 60)]


Unpack tuples:

.. code-block:: python

   bart = ('male', 10, 'simpson')  # create a tuple
   (sex, age, surname) = bart      # assign three values at once
   print(sex)
   print(age)
   print(surname)


Strings
-------

**Properties of lists: iterable, immutable**

Creating strings:

Unpack tuples:

.. code-block:: python

   s = str(42)
   s = str(1.2)
   s = str(True)
   s = 'Hello World!'

String slicing is like list slicing:


.. code-block:: python

   s[:6] # returns 'Hello '
   s[7:] # returns 'orld!'
   s[-1] # returns '!'

String methods. A string is immutable. It can not be changed in place (just like sets):

.. code-block:: python

   s.lower()           # returns 'hello world!'
   s.upper()           # returns 'HELLO WORLD!'
   s.startswith('H')   # returns True
   s.endswith('orld!') # returns True
   s.isdigit()         # returns False (returns True if every character in the string is a digit)
   str.isdigit('1.23') # Is also False, because '.' is not a digit, although that string could be a float.
   s.find('World')     # returns index of first occurrence (6), but doesn't support regex
   s.find('Planet')    # returns -1 since not found
   s.replace('Hello', 'Goodbye') # returns a string where all instances of 'Hello' are replaced with 'Goodbye'

A string can be split into lists. Default delimiter is space (' '). But different delimiters can be passed. This functionality can come in useful when parsing text files:

.. code-block:: python

   s = 'I like you'
   s.split(' ')        # returns ['I', 'like', 'you']
   s.split()           # equivalent (since space is the default delimiter)
   s2 = 'a, an, the'
   s2.split(',')       # returns ['a', ' an', ' the']

A list of strings can be joined into a single string:

.. code-block:: python

   stooges = ['larry', 'curly', 'moe']
   ' and '.join(stooges)   # returns 'larry and curly and moe'

Arithmetic operations on strings. Some arithmetic operations work on strings. Addition, for example, concatenates strings:

.. code-block:: python

   s3 = 'The meaning of life is'
   s4 = '42'
   s3 + ' ' + s4       # returns 'The meaning of life is 42'
   s5 = 'kartoffe' + 'lauf' * 2
   print(s5)

Removing whitespaces and unwanted leading/trailing characters:

.. code-block:: python

   s6 = '   spam    '
   s6.strip()         # returns 'spam'
   s7 = '$ ls && pwd'
   s7.lstrip('$ ')    # returns 'ls && pwd'


String substitutions:

.. code-block:: python

   'raining %s and %s' % ('cats', 'dogs')                       # old way
   'raining {} and {}'.format('cats', 'dogs')                   # newer way
   'raining {arg1} and {arg2}'.format(arg1='cats', arg2='dogs') # named arguments
   arg1 = 'cats'
   arg2 = 'dogs'
   f'raining {arg1} and {arg2}'                                 # newest way


String formatting:

.. code-block:: python

   # whitespace alignment
   animals = ['cat', 'dog', 'horse', 'crocodile', 'pidgeon', 'tortoise']
   for a in animals:
       print(f'The {a:<10} is an animal')

   numbers = [1, 200, -3, -5.234, 3.14159, -3.14159]
   for n in numbers:
       print(f'{n:-.2f}')
   for n in numbers:
       print(f'{n:5.3f}')
   for n in numbers:
       print(f'{n:07.2f')

Normal strings vs raw strings:

.. code-block:: python

   print('first line\nsecond line')    # normal strings allow for escaped characters
   print(r'first line\nfirst line')    # raw strings treat backslashes as literal characters

Dictionaries
------------

**Properties of unordered, iterable, mutable**

**made of key-value pairs**

**keys can be str, numbers, tuples**

**values can be anything**

Create empty dictionaries:

.. code-block:: python

   empty_dict = {}
   empty_dict = dict()
   bool(empty_dict) # is False

create a dictionary (two ways):

.. code-block:: python

   family = {'dad':'homer', 'mom':'marge', 'size':6}
   family = dict(dad='homer', mom='marge', size=6)

Convert a list of tuples into a dict:

.. code-block:: python

   list_of_tuples = [('dad', 'homer'), ('mom', 'marge'), ('size', 6)]
   family = dict(list_of_tuples)

Convert two lists into a dict:

.. code-block:: python

   keys = ['dad', 'mom', 'size']
   values = ['homer', 'marge', 6]

   family = dict(zip(keys, values))

Examing a dict:

.. code-block:: python

   family['dad']       # returns 'homer'
   len(family)         # returns 3
   'mom' in family     # returns True
   'marge' in family   # returns False (only checks keys)

Accessing keys, value and items:

.. code-block:: python

   list(family.keys())       # keys: ['dad', 'mom', 'size']
   list(family.values())     # values: ['homer', 'marge', 6]
   list(family.items())      # key-value pairs: [('dad', 'homer'), ('mom', 'marge'), ('size', 6)]

Why the extra ``list()``?

Calling ``family.keys()`` and the other dict methods doesn't return a list. It returns what is called a *view object*.

Dictionaries are mutable, so they can be changed in-place.

.. code-block:: python

   family['cat'] = 'snowball'              # add a new entry
   family['cat'] = 'snowball ii'           # edit an existing entry
   del family['cat']                       # delete an entry
   family['kids'] = ['bart', 'lisa']       # dictionary value can be a list
   family.pop('dad')                       # remove an entry and return the value ('homer')
   family.update({'baby':'maggie', 'grandpa':'abe'})   # add multiple entries

Accessing values via indexing and using the ``get()`` method:

.. code-block:: python

   family['mom']                       # returns 'marge'
   family.get('mom')                   # equivalent
   family['grandma']                   # throws an error since the key does not exist
   family.get('grandma')               # returns None instead
   family.get('grandma', 'not found')  # returns 'not found' (the default)

Accessing lists in dictionaries:

.. code-block:: python

   family['kids'][0]                   # returns 'bart'
   family['kids'].remove('lisa')       # removes 'lisa'

String substitution usind dicts:

.. code-block:: python

   'youngest child is %(baby)s' % family   # returns 'youngest child is maggie'

Sets
----

**properties: unordered, iterable, mutable, can contain multiple data types**

Sets are like dictionaries with just keys. They also use the curly braces in their definition. That's why you can't create an empty set with ``{}``, which would return an emtpy dict.

.. code-block:: python

   empty_set = set()
   bool(empty_set)   # is False

Sets with values in them:

.. code-block:: python

   languages = {'python', 'r', 'java'}         # create a set directly
   snakes = set(['cobra', 'viper', 'python'])  # create a set from a list

Examine a set:

.. code-block:: python

   len(languages)              # returns 3
   'python' in languages       # returns True

Set oprations:

.. code-block:: python

   languages & snakes          # returns intersection: {'python'}
   languages | snakes          # returns union: {'cobra', 'r', 'java', 'viper', 'python'}
   languages - snakes          # returns set difference: {'r', 'java'}
   snakes - languages          # returns set difference: {'cobra', 'viper'}

Modify a set. Sets are mutable and thus can be altered in-place:

.. code-block:: python

   languages.add('sql')        # add a new element
   languages.add('r')          # try to add an existing element (ignored, no error)
   languages.remove('java')    # remove an element
   languages.remove('c')       # try to remove a non-existing element (throws an error)
   languages.discard('c')      # remove an element if present, but ignored otherwise
   languages.pop()             # remove and return an arbitrary element
   languages.clear()           # remove all elements
   languages.update(['go', 'spark'])  # add multiple elements (can also pass a set)

Get a sorted list of unique elements from a set:

.. code-block:: python

   sorted(set([9, 0, 2, 1, 0]))    # returns [0, 1, 2, 9]

**Beware**: Sometimes the ``list(set(list(...)))`` construction is used to remove duplcate entries from lists. However, because sets are inherently unordered the resulting list might have the wrong order of elements.


Defining functions
------------------

Define a function with ``def``:

.. code-block:: python

   def print_text():
       print('this is text')

Calling the function. Calls are always executed, when an object is followed by parentheses. That's why ``a = 3; a()`` will raise the Error: ``TypeError: int is not callable``. But a function can be called with:

.. code-block:: python

   print_text()

Create a function with one positional argument:


.. code-block:: python

   def print_this(x):
       print(x)

Calling this function. If a function has a ``return`` statement, it can return objects. However, the ``print_this()`` function has no ``return`` statement and, thus, returns the default ``None``.

.. code-block:: python

   print_this(3)       # prints 3
   n = print_this(4)   # prints 4 and assigns None to n
   print(n)            # prints None

Return statements are written like this:

.. code-block:: python

   def square_this(x):
       return x**2

Adding documentation to your function. Some people might want to use your function but don't want to read the whole code of the function. You can help them with a so-called docstring (three consecutive doulbe quotes or single quotes ``"""`` or ``'''``) and give a summray of a functions code. That way it's easier to reuse the function.

.. code-block:: python

   def square_this(x):
       """Returns the square of the provided argument."""
       return x**2

Call the function or assign the return value to a new variable:

.. code-block:: python

   square_this(3)          # prints 9
   var = square_this(4)    # assigns 16 to var, but does not print 16

Positional arguments and Keyword arguments. Keyword arguments have default values and are defined like so:

.. code-block:: python

   def calc(a, b, op='add'):
       if op == 'add':
           return a + b
       elif op == 'sub':
           return a - b
       else:
           print('valid operations are add and sub')

Calling this function:

.. code-block:: python

   calc(10, 4, op='add')   # returns 14
   calc(10, 4, 'add')      # also returns 14: unnamed arguments are inferred by position
   calc(10, 4)             # also returns 14: default for 'op' is 'add'
   calc(10, 4, 'sub')      # returns 6
   calc(10, 4, 'div')      # prints 'valid operations are add and sub'

Placeholder functions. If you know you want to implement some functions, you can use ``pass`` as a placeholder (a function with an empty body is not allowed):

.. code-block:: python

   def stub():
       pass

Multiple vlues are returned as tuples:

.. code-block:: python

   def min_max(nums):
      return min(nums), max(nums)

Call this function:

.. code-block:: python

   nums = [1, 2, 3]
   min_max_num = min_max(nums)      # min_max_num = (1, 3)
   print(type(min_max_num))         # tuple
   min_num, max_num = min_max(nums) # direct unpacking

Packing and Unpacking
---------------------

Use ``*`` to unpack lists and other iterables (like tuples) into multiple arguments:

.. code-block::

   def addition_of_four(a, b, c, d):
       return a + b + c + d

   mylist = [1, 2, 3, 4]

   added_numbers = addition_of_four(*mylist) # works
   added_numbers = addition_of_four(mylist)  # won't work.

Instead of fixing the numbers of arguments to a function, you can let the function take an arbitrary amount of arguments (this is called packing). The ``*args`` is a convention. You don't need to call it ``*args``, you can also call it ``*numbers``.

.. code-block:: python

  def addition(*args):
      return sum(args)

  addition(1, 5, 10, 20)
  additon(1, -1)

That's it for unnamed (positional) arguments. But what about named arguments? For that we have ``**``:

.. code-block:: python

   def subtraction(minuend=0, subtrahend=0):
       return minuend - subtrahend

   mydict = {'minuend': 10, 'subtrahend': 5}

   subrtaction(**mydict)

Of course, when ``mydict`` has more or different keys, the function ceases to work:

.. code-block:: python

   mydict.update({'another_key': 'Hello World!'})
   subtraction(**mydict)

That's why we can also pack an arbitrary number of keyword arguments. The convention here is to use ``**kwargs``:

.. code-block:: python

   def subtraction(**kwargs):
       return kwargs['minuend'] - kwargs['subtrahend']
   subtraction(**mydict)

Lambda functions
----------------

Used to temporarily define small throwaway-functions.

.. code-block:: python

   def squared(x):
       return x**2

   squared = lambda x: x**2

Instead of this:

.. code-block:: python

   simpsons = ['homer', 'marge', 'bart']
   def last_letter(word):
       return word[-1]
   simpsons_sorted = sorted(simpsons, key=last_letter)

You can use ``lambda`` to abbreviate this to:

.. code-block:: python

   simpsons = ['homer', 'marge', 'bart']
   simpsons_sorted = sorted(simpsons, key=lambda word: word[-1])

Comprehensions
--------------

Comprehensions, like the lambda functions, are ways to make code cleaner and oftentimes easier to read. Consider this example, where a list of cubes is created:

.. code-block:: python

   nums = [1, 2, 3, 4, 5]
   cubes = []
   for num in nums:
       cubes.append(num**3)

Instead, you can also do:

.. code-block:: python

   cubes = [num**3 for num in nums]

Comprehensions can contain ``if``-``else``. Instead of:

.. code-block:: python

   cubes_of_even = []
   for num in nums:
       if num % 2 == 0:
           cubes_of_even.append(num**3)

You could also do:

.. code-block:: python

   cubes_of_even = [num**3 for num in nums if num % 2 == 0]

You can do similar things with dictionaries:

.. code-block:: python

   fruits = ['apple', 'banana', 'cherry']
   fruit_lengths = {fruit:len(fruit) for fruit in fruits}              # {'apple': 5, 'banana': 6, 'cherry': 6}
   fruit_indices = {fruit:index for index, fruit in enumerate(fruits)} # {'apple': 0, 'banana': 1, 'cherry': 2}

Map and Filter
--------------

``map`` applies a function to every element in an sequence (list, tuple, set), whereas ``filter`` uses a function that returns ``True`` or ``False`` to remove the ``False`` elements from a sequence:

.. code-block:: python

   simpsons = ['homer', 'marge', 'bart']
   list(map(len, simpsons))                      # returns [5, 5, 4]
   list(map(lambda word: word[-1], simpsons))    # returns ['r', 'e', 't']

   # equivalent list comprehensions
   [len(word) for word in simpsons]
   [word[-1] for word in simpsons]

Again, the ``list()`` is because ``map`` returns a ``iterator`` which can not directly be printed. (Why does python do this? It postpones the actual execution of the mapping function until it is needed. This saves computation time if you don't need the ``map`` applied to the whole list but, for example, use it in a for loop, that stops early). Finally, here's ``filter``:

.. code-block:: python

   nums = range(5)
   filter(lambda x: x % 2 == 0, nums)      # returns [0, 2, 4]

   # equivalent list comprehension
   [num for num in nums if num % 2 == 0]

Static quick reference notebooks
================================

.. toctree::
   :maxdepth: 1

   quick_references/qr_01_data_types.nblink

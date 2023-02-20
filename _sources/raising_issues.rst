 .. _raising-issues-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

==============
Raising Issues
==============

The *raise issue* functionality on GitHub can be used for comments, remarks or questions. If you find yourself stuck at a problem and it won't budge head over to GitHub and raise an issue. Here, I will show you how to raise an issue, that references a specific line in the python_basics introduction. From the main github page (https://github.com/kevinsawade/start-science-here) I make my way to the file containing the line I want to reference in a new issue. I go to ``python_tutorial/01_python_basics.ipynb``. If you can't find a file inside a GitHub repository, you can try the search bar in the top left. However, I know where the file is and want to raise an issue from a line in that file. This file is special, because GitHub renders the ipython notebook file (.ipynb), just like it would render .pdf, .png and r.st files (Have a look at https://github.com/kevinsawade/start-science-here/blob/main/docs/source/index.rst to see a rendered .rst file). I will enter the unrendered raw view via the <> button:

.. image:: _static/pics/raising_and_merging/raise_issue_1.PNG
   :target: _static/pics/raising_and_merging/raise_issue_1.PNG
   :alt: Image showing how to view a raw file.

Once in the raw view, I will raise an issue from a line by clicking the three dots at the start and selecting "reference in new issue".

.. image:: _static/pics/raising_and_merging/raise_issue_2.PNG
   :target: _static/pics/raising_and_merging/raise_issue_2.PNG
   :alt: Image showing how to raise issue from line.

I will then give the issue a descriptive name and some text. This window is basically the same, as if you were raising a simple issue, but it references the line I raised the issue from via a link.

.. image:: _static/pics/raising_and_merging/raise_issue_3.PNG
   :target: _static/pics/raising_and_merging/raise_issue_3.PNG
   :alt: Image showing how to compose an issue.

After that I will wait for someone to address my issue. In this case, I will answer my own issue from another account. :)

.. image:: _static/pics/raising_and_merging/raise_issue_4.PNG
   :target: _static/pics/raising_and_merging/raise_issue_4.PNG
   :alt: How an issue is resolved.

The issue has been resolved via a new git commit that implements a solution to my issue.

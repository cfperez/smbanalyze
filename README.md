# smbanalyze

An analysis package to analyze single-molecule optical trapping and FRET data. The movitation: I like Python, and I hate IgorPro.

The original implentation was very heavy-handed--bloated classes, unnecessary classes, extensive heirarchies. That the current state of the code is an *improvement* is kind of embarrassing, but it is a journey after all!

If I were to do it all over, I would either 1) use pandas or 2) create a lighter-weight skin over numpy arrays. As for option 1, my purpose in this project (besides creating usable software) was to experience and better understand the process of creating a software suite with actual users. How do I make the APIs? Name things? Relate this object to that one? Split code into modules? I may not have learned these things if I hopped into pandas, although I might have gained exposure to *good* design upfront--instead I was influenced by old C++ knowledge and matplotlib.

What I noticed at the end is that our data are essentially just numpy arrays with metadata. Doing it all over (and not using pandas), I would spend the time to learn how to properly subclass nd.array and simply add a .metadata attribute dictionary instead of using my own custom Datatype. However, I did learn some neat Python metaprogramming techniques!

## Usage

Load the package for interactive use:
```python
    from matplotlib.pyplot import *
    from smbanalyze.shell import *
```
Gives you immediate access to most packages plus some convenience functions. This is the recommended way of using the package for interactive use.

See notebooks/ for (somewhat hacky) example code.

experiment.py and datatypes.py define the data types experiment.Pulling, TrapData, and FretData (these last two may by simplified in the future.)

fec.py fits WLC to Pulling experiments and extracts features.

refolding.py is used for calculating binding curves from refolding time experiments.

db/ is used for storing/retrieving experiments saved in databases. Currently, this only houses mongoDB functions, but in the future, it can be used to create a local SQLlite database for potentially more robust data crunching (if needed.)

fplot.py has the plotting functions for pretty graphs.

fcalc.py calculates the fret from the image data.

### TODO/needs

+ Docstrings! Duh.
+ Unit tests -- especially for regression. It is currently used for real analysis, so I know its working!
+ Eliminate cruft -- fplot, fileIO, and numlist are probably the worst offenders.
+ FRET analysis -- Bayesian state assignment, etc.
+ DB analysis -- complex analysis on stored data
+ Better interactive graphics -- matplotlib is difficult to use to create truly interactive applications
+ Save session -- automatically save certain objects or classes to a local persistent dictionary or DB (per session)
+ Project organization -- update to best practices for python packages, including tox, CI, distribute, etc.

### What I've learned

#### Functional over object oriented

While I'm not doing functional programming, my programming style has become much more functional. I now understand better  how (and why) software is written to create loose coupling or orthoganol components.

If new behaviors, functions, or processes are often needed but the data stays the same, then program functionally.
If new objects or types are often needed but the behaviors are the same, then program "objectively".

I consider the program and problem needs more abstractly now. This is especially relevant to data science where the data type system is fairly standard--strings, numbers, lists, hashes, graphs--and its transformations, operations, and analyses of these data that are constantly in flux.

I'm still learning when to subclass or compose objects instead, but I will pursue a functional solution first unless/until the program complexity/maintainability would benefit from a new type.

#### Get users

I cared a lot because I was using this to do actual work. I cared even more because my labmates were using it. I wanted to fix errors ASAP, and learn how to make better code so I could avoid them in the first place. Otherwise, the experience can be very academic and I fall into the trap of making some unimportant function really "perfect".

#### Use the minimum viable objects

Create a data type that only exposes methods that are fundamental to the behavior of the object [a la Code Complete]. Often, this can be done with Python's excellent built-in types or stdlib collections. Using them means objects are guaranteed to work and fit together with all the other parts of the language.

I struggle (struggled?) with making things too OOP-y and brittle to changes. While frustruating, my search for answers brought me to Google talks on The Law of Demeter, composition, and testability. So when objects are required, I should be able to make them better (more testable, less brittle) than before.

#### API design is hard

Especially without TDD and regression tests to help you go backwards, but even then! In the end, I think I finally came to an understanding about how to split up the code to reflect the different conceptual components of the problem domain. Understanding the problem space helps create code that is modularized in a way that makes sense to the humans using it.

#### Learn the fundamentals

I will tackle fundamentals in CS in the same way and for the same reason that I approached physics--for general principles that are widely applicable for creating things I care about. For me, I think of CS as all about methods to handle complexity in information through abstraction and math. Unsurprisingly, this is how I think of physics--how to handle the complexity in nature through abstraction and math. Because I'm expert at learning principles and applying them to problems, I know that studying CS fundamentals will be much more than a requirement to passing interviews. I will use it to improve my ability to create things in software.

#### How to work without saved state?

In IPython Notebook, if the instance crashes, all the local variables are gone. This is bad when you are doing manual analysis of single-molecule traces and are in the middle of some procedure. Unlike other environments, I found no easy/transparent way to save state of my current progress. Maybe I'm ignorant of some obvious solution, or I am mistakenly assuming that I should be able to work that way!

#### Role models and code review

I wish I had more of both. Programming is a craft built on fundamental techniques. Techniques can be learned in a book and practiced in isolation. After that, good judgement--when, how, and why to use them--becomes critical, and good judgement comes from experience. Working with a superior programmer would allow me to learn from their experience by closely examining how they make choices.

#### TDD?

Definitely open to trying this again. I gave it a shot, but my code was so badly formed at my first early attempts that I was getting mired in mocks and patches. Later attempts were easier because the code was better designed... and worked so much better that I felt ok not writing tests unless I wanted to check a particular corner case.

#### Github

Obviously :)

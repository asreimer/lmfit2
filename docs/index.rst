lmfit2
======

**lmfit2** is an GPLv3 licensed implementation of the First-Principles
Fitting Methodology (FPFM) algorithm described in `Statistically Self‐Consistent
and Accurate Errors for SuperDARN Data<https://doi.org/10.1002/2017RS006450>`_.
This algorithm was specifically developed for fitting SuperDARN "rawacf" data,
however the principles used to develop the algorithm are broadly applicable to
parameter extraction via fitting models to data. lmfit2 was written in two
different languages, C and Python, but the C version is the recommended version
to use. These pages will show you how to use lmfit2.

This documentation will not teach you much about how to fit data, but there are
many good resources on this topic already available (try `numerical recipes
<http://numerical.recipes/>`_ specifically the chapter on modeling data).
We also `published a paper <https://doi.org/10.1002/2017RS006450>`_ explaining
the FPFM algorithm and implementation in detail.

lmfit2 fitted data was used `in a study identifying ULF waves in SuperDARN data
<https://doi.org/10.1029/2018JA025859>` and it is being developed on `GitHub
<https://github.com/asreimer/lmfit2>`_.

.. image:: https://img.shields.io/badge/GitHub-asreimer%2Flmfit2-blue.svg?style=flat
    :target: https://github.com/asreimer/lmfit2
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg?style=flat
    :target: https://github.com/asreimer/lmfit2/blob/master/LICENSE



Basic Usage
-----------

For the C version, if you want to fit a rawacf file (in dmap format), you would
do something like this:

.. code-block:: bash

    make_lmfit2 -new yyyymmdd.hhmm.rawacf > yyyymmdd.hhmm.lmfit2


For the python version, with `lmfit2.py` in the working directory, you would do
something like this:

.. code-block:: python

    from lmfit2 import main
    main(yyyymmdd.hhmm.rawacf,yyyymmdd.hhmm.lmfit2)


How to Use This Guide
---------------------

First, you need to install lmfit2. Luckily for you there's a handy
:ref:`install` guide to help you.

Next, you can try using lmfit2 as explained above in the Basic Usage section.

Perhaps some tutorials will be added to this guide in the future! If you would
like this to happen, or you have bug reports, patches, feature requests, and/or
other comments, please submit them to `the GitHub issue tracker <https://github.com/dfm/emcee/issues>`_.

If you have a question about the use of lmfit2, please post it to the issue tracker.


.. toctree::
   install
   ussage


License & Attribution
---------------------

Copyright 2016-2020 Ashton S. Reimer.

lmfit2 is free software made available under the GPLv3 License. For details
see the ``LICENSE``.

If you make use of lmfit2 in your work, please cite our paper
(`Radio Science <https://doi.org/10.1002/2016RS005975>`_,
`ADS <https://ui.adsabs.harvard.edu/abs/2016RaSc...51..690R/abstract>`_,
`BibTeX <https://ui.adsabs.harvard.edu/abs/2016RaSc...51..690R/exportcitation>`_)


Changelog
---------

.. include:: ../HISTORY.rst

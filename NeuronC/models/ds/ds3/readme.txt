
We have a DS model!  Let's share it -

Here's how:  Everyone checks out the model at their own pace --
when you need something modified, you request modifications.
I'll do my best to modify the script per your requests -- but of
course I'd like to encourage you to work on the script yourself.
If you request a private modification, I'll send the modification
back to you only, but otherwise I'll let everyone on the DS list
have it.  The idea is that if everyone participates, then we all
will gain from it.  There is a learning curve for the simulator
and the model script, but I believe that the time that we spend
on learning will pay back in understanding. 

To run the model, you'll need to download the simulator and
compile it -- for this, you'll need a Linux (or any Unix)
machine.  There is also a DOS version but I don't recommend it,
because the Linux environment is so much more productive.  A
manual in .html form is included in the simulator distribution.
You can get it from:

    ftp://retina.anatomy.upenn.edu/pub/nc.tgz

The DS model is attached as a .zip file (the package contains
several scripts and some auxiliary files).  The basic script
("ds3.n") describes the morphology, connections, and experiments,
and there is also a small library of functions ("gcdefs.n"). The
model has one ganglion cell and an array of bipolars and
amacrines feeding it with synaptic inputs.  There is no feedback,
only feedforward, and the inputs are very simple -- a direct
translation of light intensity into a voltage in the bipolar
cell.

The model has only "generic" morphology for the GC and amacrines
which is constructed by an algorithm.  The bipolars are just
single compartments.  But there is a provision for adding
realistic morphology for either cell type.  If someone would like
to provide cell morphology for the public domain I'll include it.
Scripts are available to convert from Neurolucida format and it's
fairly easy to convert from other formats.  The format is defined
in the manual and in the "gcdefs.n" script.

An array of bipolars receives signals from a moving bar stimulus.
The bipolars are randomly located, specified by a
nearest-neighbor distance mean (bp_nn) and regularity
(mean/stdev. = 1/coefficient of variation, default 10).  The
bipolars provide synaptic input to the ganglion cell and to an
array of amacrines.  The bp -> gc synapses are created whenever a
bipolar cell is within a criterion distance of a ganglion cell
dendrite ("bpsyn_dist_thresh", set by default to 7 um).

The amacrine cells are also located randomly with a
nearest-neighbor distance mean (am_nn) and regularity (def. 5).
The default amacrine morphology is a soma and 1 dendrite,
although the algorithm has a provision to include multiple
dendrites.  A bipolar cell provides synaptic input to an amacrine
when the 2 cells are within a criterion distance.

An amacrine provides synaptic input to the ganglion cell when
cable segments of the 2 cells are within a criterion distance,
and the synapse would be farther than a criterion distance away
from previous synapses on the same amacrine, and also a criterion
distance from the amacrine soma.

I hope that this model will generate some new ideas about
possible mechanisms -- because to add new mechanisms we will need
to communicate -- and I expect that this communication will be
helpful to everyone.

Please let me know if you have comments or requests for changes
on the script -- or advice about how we should proceed or
communicate.   

Best wishes,

Rob

rob@retina.anatomy.upenn.edu
Aug, 2001


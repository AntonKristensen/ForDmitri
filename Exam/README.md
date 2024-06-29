Hi Dmitri!

I have number 32 on the list, so my exam task was:
"Symmetric rank-1 update of a size-n symmetric eigenvalue problem"

I have solved it, mostly in the exam.cpp file. You can see the result of the timing on ./plotting/time.png.
I have fitted an a*x^2 function to it, and it matches decently.

It sometimes has trouble finding the roots. If I allow the root finding algorithm (in root.cpp) to run for a long time, these outliers will give lots of extra timing, as can be seen on "TimeWithExaggeratedOutliers.png". I have currently put it to 50. I tested and didn't see any normal runs taking more than 10 counts.

I did manage to fix most of the occasions when it couldn't find the root. It turned out that the root was close to d[i], so it needs a precise guess
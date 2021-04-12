Daqana VERSION 4 "daqana_v6f4" was implemented to make a pulse processing
with an average in the the baseline before the rising time of the pulse and in the end of it
and not around the minimum.

The reason why is implemented like this, is because we are having problems 
with an overshot in the minimum of the pulse. We at the moment don't
understand what is it causing this but we would try to skip this by 
calculating the average a little further away the minimum where this over
shut has came back to the real value of the post baseline already.

201216, Ana.

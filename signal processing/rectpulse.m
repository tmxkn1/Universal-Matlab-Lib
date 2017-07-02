function f = rectpulse(t,t1,t2)
f = (t-t1>=0) - (t-t2>=0);
end
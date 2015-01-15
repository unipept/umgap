rm _rmq.so rmq.o rmq.py rmq.pyc rmq_wrap.c rmq_wrap.o
swig -python rmq.i
cc -c `python3-config --cflags` -fPIC rmq.c rmq_wrap.c
ld -shared rmq.o rmq_wrap.o -o _rmq.so

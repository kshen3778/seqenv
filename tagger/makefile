CC     = g++
CFLAGS = -fpic -Wall -O3
LFLAGS = -fpic -shared -lboost_regex -lpython2.7
PYTHON = -I/usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/include/python2.7 -I/usr/local/Cellar/boost/1.60.0_1/include/boost/tr1

all: tagger_swig.py _tagger_swig.so libtagger.so

clean:
	rm -f tagger_swig.py *.pyc *_wrap.cxx *.o *.so

#
# Make C-style Python interpreter wrapper C-code layer.
#
tagger_swig.py tagger_swig_wrap.cxx: tagger_swig.i tagger.h
	swig -python -c++ -threads $<

tagger.o: tagger.cxx acronyms.h tagger.h tagger_types.h tagger_core.h acronyms.h hash.h tokens.h
	$(CC) $(CFLAGS) -c $<

tagger_swig_wrap.o: tagger_swig_wrap.cxx tagger.h tagger_types.h tagger_core.h acronyms.h hash.h tokens.h
	$(CC) $(CFLAGS) $(PYTHON) -c $<

#
# Link to make the Python wrapped C/C++ tagger interface in
# the Shared Object: _tagger_swig.so.
#
_tagger_swig.so: tagger_swig_wrap.o
	$(CC) $(LFLAGS) -o $@ $<

#
# Link to make the C/C++ stand-alone Shared Object for
# inclusing the C/C++ tagger in other C++ programs.
#
libtagger.so: tagger.o
	$(CC) $(LFLAGS) -o $@ $<

%: %.cxx acronyms.h document.h file.h hash.h mutex.h match_handlers.h base_handlers.h batch_tagger.h tagger.h tagger_core.h tagger_types.h tightvector.h tokens.h
	$(CC) $(CFLAGS) -lboost_regex -o $@ $< -lm

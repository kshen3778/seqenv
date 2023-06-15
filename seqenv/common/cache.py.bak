# Built-in modules #
import inspect

################################################################################
def cached(f):
    """Decorator for functions evaluated only once."""
    def memoized(*args, **kwargs):
        if hasattr(memoized, '__cache__'):
            return memoized.__cache__
        result = f(*args, **kwargs)
        memoized.__cache__ = result
        return result
    return memoized

################################################################################
def property_cached(f):
    """Decorator for properties evaluated only once.
    It can be used to created a cached property like this:

        class Employee(object):
            @property_cached
            def salary(self):
                print "Evaluating..."
                return time.time()
        bob = Employee()
        print bob.salary
        bob.salary = "10000$"
        print bob.salary
    """
    def retrieve_from_cache(self):
        if '__cache__' not in self.__dict__: self.__cache__ = {}
        if f.__name__ not in self.__cache__:
            # Add result to the property cache #
            if inspect.isgeneratorfunction(f): result = tuple(f(self))
            else: result = f(self)
            self.__cache__[f.__name__] = result
        return self.__cache__[f.__name__]
    def overwrite_cache(self, value):
        if '__cache__' not in self.__dict__: self.__cache__ = {}
        self.__cache__[f.__name__] = value
    retrieve_from_cache.__doc__ = f.__doc__
    return property(retrieve_from_cache, overwrite_cache)

###############################################################################
class class_property(property):
    """You can use this like this:
        class Foo(object):
            @class_property
            @classmethod
            @cached
            def bar(cls):
                return 1+1
    """
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()
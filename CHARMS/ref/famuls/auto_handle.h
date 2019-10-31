#ifndef AUTO_HANDLE_H
# define AUTO_HANDLE_H

// A template class for returning auto_handles by value
template <class H>
struct auto_handle_ref
{
  explicit auto_handle_ref (H h) : _handle(h) {}
  H _handle;
};

// auto_handle
template<class H,
         class T = handle_traits<H>,
         class B = base_handle<H,T> >
class auto_handle : private B {
    public:
    // construct/copy/destroy
    explicit auto_handle (const H& h = T::null_value()) throw() : _handle(h) {}

    auto_handle(auto_handle &that) throw () : B(that), _handle(B::copy(that._handle)) {}
    ~auto_handle() throw ()
    {
      B::dispose (_handle);
    }

    // members
    H get() const throw() { return _handle; }

    // release ownership
    H release () throw()
    {
      H h(_handle);
      _handle = T::null_value();
      return h;
    }

    void reset (const H &handle) throw()
    {
      if (_handle != handle) {
        T::dispose (_handle);
        _handle = handle;
      }
    }

    // conversions
    // implicit ctor, clients may write auto_handle<some_class> h = func()
    // where func returns auto_handle by value
    auto_handle(const auto_handle_ref<H> &r) throw() : _handle(r._handle) {}

    operator auto_handle_ref<H>() { return auto_handle_ref<H>(release()); }

    // other operators
    auto_handle &operator=(const auto_handle_ref<H> &r)
    {
      auto_handle tmp(r);
      std::swap(_handle, tmp._handle);
      return *this;
    }

    auto_handle &operator=(auto_handle_ref<H> &r)
    {
      auto_handle tmp(r);
      std::swap(_handle, tmp._handle);
      return *this;
    }

    bool operator!() const { return _handle == T::null_value(); }

    private:
    H _handle;
  };

#endif

////////////////////////////////////////////////////////////////
// $Id: Cube::Handle.hxx,v 1.34 2013/01/11 16:52:32 mcgrew Exp $
//
#ifndef CubeHandle_hxx_seen
#define CubeHandle_hxx_seen

#include <iostream>
#include <typeinfo>

#include <TObject.h>

#include "CubeLog.hxx"

/// An implementation of a shared_ptr (like) object that can be
/// streamed to a TFile by ROOT.
namespace Cube {

    class VHandle;
    class HandleBase;
    class HandleBaseDeletable;
    class HandleBaseUndeletable;

    template <class T> class Handle;
    template <class T> T* GetPointer(const Handle<T>& handle);

    template <class T, class U> bool operator < (const Handle<T>& a,
                                                 const Handle<U>&b);

    /// An abstract base class for handles that's used to maintain the
    /// reference count.  This is a deep internal class that can't be directly
    /// accessed.  The original didn't define a constructor, or destructor so
    /// it won't interfere with templated objects that derive from this class,
    /// but they are needed by root so this class has them...
    class VHandle: public TObject {
    protected:
        VHandle();
        virtual ~VHandle();

        /// Define default values for this object.
        void Default(HandleBase* handle);

        /// Add a reference to the object being held by adding this
        /// Cube::Handle to the reference list.
        void Link(const VHandle& rhs);

        /// Remove a reference to the object being held by removing this
        /// Cube::Handle from the reference list.  Returns true if the last
        /// reference is removed.
        void Unlink();

        /// Safely get the pointer value for this handle.  This hides the
        /// underlying storage model from the Cube::Handle template.
        TObject* GetPointerValue() const;

        /// Define the status bits used by the VHandle object.  These can't
        /// collide with any status bits defined in TObject (the parent class
        /// for VHandle), and none of the VHandle children can define a
        /// status bit that collides with these definitions (i.e. the
        /// Cube::Handle<> templates.  Bits 14 to 23 are available for use.
        enum EStatusBits {
            kWeakHandle = BIT(20)
        };

    public:
        /// Release the ownership of the object being held by this handle.
        /// The responsiblity to delete the object passes to the calling
        /// routine.
        void Release();

        /// Make the current handle into a weak handle.  A weak handle does
        /// not own the object and the object may be deleted.  If the object
        /// is deleted, the weak handle will reference a NULL pointer.  The
        /// MakeWeak method may remove the last reference to the object and
        /// cause it to be deleted, so it should be check for validity before
        /// using.
        void MakeWeak();

        /// Make the current handle into a regular handle that "owns" the
        /// object.  The object won't be deleted until all handles that own
        /// the object are removed.
        void MakeLock();

        /// Check if this is a weak pointer to the object.
        bool IsWeak() const {return TestBit(kWeakHandle);}

        /// Equality operator for all Cube::Handle objects.
        bool operator == (const VHandle& rhs) const;

        /// A deep debugging tool, DO NOT CALL!!!!
        HandleBase* GetInternalHandle() const {return fHandle;}

        /// Print the hit information.
        virtual void ls(Option_t *opt = "") const;

    private:
        /// Check if the object, or the handle to the object should be
        /// deleted.  This can be safely called even if the object and handle
        /// don't exist.
        void CheckSurvival();

        /// The reference counted handle. This handle contains the pointer to
        /// the actual data object.
        HandleBase* fHandle;

        ClassDef(VHandle,10);
    };

    /// A handle class that will manage the lifetime (and ownership) of
    /// objects that are derived from TObject.  This doesn't require any
    /// cooperation from the object that will be reference counted.  This
    /// works fine even if the object is owned by ROOT, but really shines if
    /// the object is owned by the program.  The resulting handle is a smart
    /// pointer and (in most respects) looks like a normal pointer.
    ///
    /// \warning Cube::Handle objects which are used in an event loop \b must \b
    /// not be saved between events.  Saving a Cube::Handle between events will
    /// cause a memory leak when the output file is read by a later program.
    /// In addition, keep in mind that this is a reference counted object (see
    /// Wikipedia for details).  This means that "reference loops" will
    /// prevent objects from being deleted.
    ///
    /// Reference loops can be prevented by using weak handles.  A weak handle
    /// is a handle that can refer to the object without increasing the object
    /// reference count.  A handle is marked as weak using the MakeWeak()
    /// method, and will remain weak until a call to MakeLock().  Note: The
    /// "weak" property follows the handle, not the object being referenced by
    /// the handle.
    ///
    /// \code
    /// Cube::Handle<Cube::Hit> aHandle(aHitFromSomePlace);
    /// aHandle.IsWeak()     // Returns false.
    ///
    /// Cube::Handle<Cube::Hit> bHandle = aHandle;
    /// bHandle.IsWeak();    // Returns false.
    /// bHandle.MakeWeak();
    /// bHandle.IsWeak();    // Returns true.
    /// bHandle.MakeLock();
    /// bHandle.IsWeak();    // Returns false.
    ///
    /// Cube::Handle<Cube::Hit> weakHandle;
    /// weakHandle.MakeWeak();
    /// weakHandle.IsWeak(); // Returns true;
    /// weakHandle = aHandle;
    /// weakHandle.IsWeak(); // Returns true;
    ///
    /// aHandle.IsWeak();    // Returns false;
    /// \endcode
    ///
    /// When handling weak and non-weak handles, the difference between
    /// assignment and the copy constructor must be observed.  The assignment
    /// operator assigns the object referenced by a handle to another handle,
    /// while the copy constructor copies the value of a handle to a new
    /// handle.  The weak-ness of the handle follows the handle and so is not
    /// transferred during assignment.  However, when a copy constructor is
    /// used, a new handle is being created, and the "weak" property of the
    /// original handle will be transfered to the new handle.  Contrast the
    /// two situations.
    ///
    /// \code
    /// Cube::Handle<Cube::Hit> aHandle(aHitFromSomePlace);
    /// aHandle.IsWeak()     // Returns false.
    ///
    /// Cube::Handle<Cube::Hit> bHandle;
    /// bHandle.IsWeak();    // Returns false.
    /// bHandle.MakeWeak();
    /// bHandle = aHandle;   // The bHandle remains weak.
    /// bHandle.IsWeak();    // Returns true.
    ///
    /// Cube::Handle<Cube::Hit> dHandle;
    /// dHandle = bHandle;   // Assign a weak to a non-weak.
    /// dHandle.IsWeak();    // The non-weak remains non-weak.
    ///
    /// Cube::Handle<Cube::Hit> eHandle(bHandle);
    /// eHandle.IsWeak();    // Returns true.
    ///
    /// Cube::Handle<Cube::Hit> fHandle(aHandle);
    /// fHandle.IsWeak();    // Returns false.
    /// \endcode
    ///
    /// An object will be deleted as soon as the last non-weak handle goes out
    /// of scope.  This has the interesting side effect that a call to
    /// MakeWeak() might remove the last non-weak handle (by making it weak),
    /// and cause the object to be deleted.  You should always check that a
    /// weak handle is valid before dereferencing it's pointer.
    ///
    /// If a weak handle can be promoted to a regular handle using the
    /// MakeLock() method.  If the weak handle holds a valid reference to an
    /// object before the call to MakeLock(), it will be a "reference counted"
    /// owner of the object after the call to MakeLock().
    ///
    /// A reference to a NULL Cube::Handle will throw an EHandleBadReference
    /// which means it won't generate a core dump.  For debugging, you can run
    /// the problem under gdb, and set a break point for the
    /// EHandleBadReference constructor.
    ///
    /// \code
    /// catch throw 'EHandleBadReference::EHandleBadReference()'
    /// \endcode
    template <class T>
    class Handle : public VHandle {
        template <class U> friend class Handle;
        template <class U> friend U* GetPointer(const Handle<U>& handle);

    public:
        /// Allow a null handle to be constructed.
        Handle();

        /// Explicitly construct a Cube::Handle from a T pointer.  If this isn't
        /// explicit, then C++ will use this constructor as a type conversion
        /// operator.  Making it explicit gives the following behavior:
        /// \code
        /// Pointee* p = new Pointee;
        /// Cube::Handle<Pointee> h(p);     // O.K.
        /// Cube::Handle<Pointee> g = h;    // O.K.
        /// Cube::Handle<Pointee> j = p;    // won't compile
        /// \endcode
        /// While conversion from a pointer to a handle good on "paper", it's
        /// not.  Internally, "h", and "g" share the ownership of "p" and
        /// manage so that "p" remains active until both "g" and "h" are
        /// destroyed.  On the other hand, "h" and "j" both assume full
        /// ownership of "p" which gets deleted as soon as the first one gets
        /// deleted (leaving a dangling reference, a.k.a. BUG).
        explicit Handle(T* pointee);

        /// Create a Cube::Handle for an object with an explicit ownership.  This
        /// allows Cube::Handle objects to refer to static objects, objects that
        /// are allocated in a buffer, or objects allocated on the stack.  If
        /// the object is owned by the handle (and can be deleted), the owner
        /// argument is true.  If the object is *not* owned by the handle,
        /// the owner argument is false.
        Handle(T* pointee, bool owner);

        /// The copy constructor for this handle.
        Handle(const Handle<T>& rhs);

        // Copy between classes.
        template <class U> Handle(const Handle<U>& rhs);

        /// The destructor for the Cube::Handle object which may delete the
        /// pointer.
        virtual ~Handle();

        /// @{ Assign one Cube::Handle object to another.  This should be
        /// designed to recast the pointee between the assignments so that an
        /// implicit conversion takes place.  If the recast fails, the new
        /// pointer will be null.
        ///
        /// This provides both const and non-const versions of the assignment.
        Handle<T>& operator = (Handle<T>& rhs);
        const Handle<T>& operator = (const Handle<T>& rhs);
        template <class U>
        Handle<U>& operator = (Handle<U>& rhs);
        template <class U>
        const Handle<U>& operator = (const Handle<U>& rhs);
        /// @}

        /// The reference operator
        T& operator*() const;

        /// The redirection operator
        T* operator->() const;

#ifndef __CINT__
    private:
        /// Internal \internal
        /// Prevent the delete operator, while allowing
        /// \code
        /// Cube::Handle<Pointee> h
        /// if (h) return true;
        /// \endcode
        class Tester {
            void operator delete(void*);
        public:
            void bogo() {}
        };

    public:
        /// A conversion to a bogus type.  This lets code like this
        ///
        /// \code
        /// Cube::Handle<Pointee> h(new Pointee);
        /// if (h) return true;
        /// \endcode
        ///
        /// compile and work as expected.  The trade off for this behavior is
        /// that we can't implement converters to the pointer type:
        ///
        /// \code
        /// Cube::Handle<Pointee> h(new Pointee);
        /// Pointee *p = h;                  // Won't work.
        /// Pointee *p = GetPointer(h);      // Will work.
        /// \endcode
        ///
        /// since it would introduce ambiguities and the compiler won't pick a
        /// method.  Actually, that's probably a benefit, The conversion to a
        /// pointer should be done with GetPointer();
        operator Tester*() const {
            if (!GetPointerValue()) return NULL;
            static Tester test;
            return &test;
        }
#endif

        ClassDefT(Handle,1);
    };
    ClassDefT2(Handle,T)

    /// Turn a Cube::Handle object into a pointer.  This is implemented as a
    /// function so that it is *really* clear that something funky is going
    /// on.  You should always pass a Cube::Handle, or a Cube::Handle reference.
    template <class T>
    T* GetPointer(const Handle<T>& handle) {
        TObject* object = handle.GetPointerValue();
        if (!object) return NULL;
        T* pointer = dynamic_cast<T*>(object);
        return pointer;
    }

    /// Make a comparision between two handles based on the pointer value.
    template <class T, class U>
    bool operator <(const Handle<T>& a, const Handle<U>& b) {
        T* aPtr = GetPointer(a);
        U* bPtr = GetPointer(b);
        return (aPtr < bPtr);
    }


    /// An abstract base class to implement the reference counted internal
    /// object.  The Cube::HandleBase objects contain the actual pointer that
    /// is being reference counted.  When the Cube::HandleBase object is
    /// deleted, the pointer is also deleted.  This object maintains the
    /// reference count.
    class HandleBase : public TObject {
    public:
        HandleBase();
        virtual ~HandleBase();

        int GetReferenceCount() const {return fCount;}
        int GetHandleCount() const {return fHandleCount;}
        void CheckHandle() {if (fHandleCount<fCount) fHandleCount=fCount;}

        // Increment/decrement the count of objects that own the object.  This
        // doesn't include any weak references.
        void IncrementReferenceCount() {++fCount;}
        void DecrementReferenceCount() {if (fCount>0) --fCount;}

        // Increment/decrement the count of objects referencing this
        // Cube::HandleBase object.  This includes the count of handles owning
        // the object as well as the weak handles that are referencing the
        // object, but don't own it.
        void IncrementHandleCount() {
            ++fHandleCount;
            if (fHandleCount > 30000) {
                CUBE_ERROR << "To many handles for object: "
                           << fHandleCount
                           << std::endl;
            }
        }
        void DecrementHandleCount() {if (fHandleCount>0) --fHandleCount;}

        // Return the current pointer to the object.
        virtual TObject* GetObject() const = 0;
        // Delete the object.  This should check that fObject is a valid
        // pointer (e.g. not NULL) that can be deleted before freeing the
        // memory.  The fObject pointer should always be set to NULL.
        virtual void DeleteObject() = 0;
        void Release() {SetBit(kPointerReleased);}
        bool IsOwner() {return !TestBit(kPointerReleased);}

    private:
        /// Define the status bits used by the Cube::HandleBase object.  These
        /// can't collide with any status bits defined in TObject (the parent
        /// class for Cube::HandleBase), and none of the Cube::HandleBase
        /// children can define a status bit that collides with these
        /// definitions.  Bits 14 to 23 are available for use.
        enum EStatusBits {
            kPointerReleased = BIT(20)
        };

        /// The number of references to the object.
        unsigned short fCount;

        /// The number of references to the handle.
        unsigned short fHandleCount;

        ClassDef(HandleBase,3);
    };

    /// A concrete version of the Cube::HandleBase class for pointers that
    /// should be deleted when the last reference goes away.  The reference
    /// count is maintained in Cube::HandleBase.
    class HandleBaseDeletable : public HandleBase {
    public:
        HandleBaseDeletable();
        HandleBaseDeletable(TObject* object);
        virtual ~HandleBaseDeletable();

        TObject* GetObject() const {return fObject;}
        void DeleteObject();

    private:
        /// The actual pointer that will be reference counted.
        TObject* fObject;

        ClassDef(HandleBaseDeletable,2);
    };

    /// A concrete version of the Cube::HandleBase class for pointers that
    /// cannot be deleted when the last reference goes away.  The reference
    /// count is maintained in Cube::HandleBase.
    class HandleBaseUndeletable : public HandleBase {
    public:
        HandleBaseUndeletable();
        HandleBaseUndeletable(TObject* object);
        virtual ~HandleBaseUndeletable();

        TObject* GetObject() const {return fObject;}
        void DeleteObject();

    private:
        /// The actual pointer that will be reference counted.
        TObject* fObject;

        ClassDef(HandleBaseUndeletable,2);
    };

    bool CleanHandleRegistry(bool);
} //End of namespace Cube.

#ifndef __CINT__
//////////////////////////////////////////////////////////////////
// Implementation of methods.
//////////////////////////////////////////////////////////////////
template <class T>
Cube::Handle<T>::Handle(T* pointee) {
    Cube::HandleBase *base = NULL;
    if (pointee) base = new Cube::HandleBaseDeletable(pointee);
    Default(base);
}

template <class T>
Cube::Handle<T>::Handle() {
    Default(NULL);
}

template <class T>
Cube::Handle<T>::Handle(T* pointee, bool owner) {
    if (pointee) {
        if (owner)
            Default(new Cube::HandleBaseDeletable(pointee));
        else
            Default(new Cube::HandleBaseUndeletable(pointee));
    }
    else {
        Default(NULL);
    }
}

template <class T>
Cube::Handle<T>::Handle(const Cube::Handle<T>& rhs) : VHandle(rhs) {
    Default(NULL);
    Link(rhs);
    if (rhs.IsWeak()) MakeWeak();
}

template <class T>
template <class U>
Cube::Handle<T>::Handle(const Cube::Handle<U>& rhs) {
    Default(NULL);
    if (dynamic_cast<T*>(rhs.GetPointerValue())) {
        Link(rhs);
    }
    if (rhs.IsWeak()) MakeWeak();
}

template <class T>
Cube::Handle<T>::~Handle() {
    Unlink();
}

template <class T>
Cube::Handle<T>& Cube::Handle<T>::operator = (Cube::Handle<T>& rhs) {
    if (operator == (rhs)) return rhs;
    // Going to replace the value of this smart pointer, so unref and
    // possibly delete.
    Unlink();
    // Compatible types
    if (dynamic_cast<T*>(rhs.GetPointerValue())) {
        Link(rhs);
    }
    return rhs;
}

template <class T>
const Cube::Handle<T>& Cube::Handle<T>::operator = (const Cube::Handle<T>& rhs) {
    // Going to replace the value of this smart pointer, so unref and
    // possible delete.
    Unlink();
    // Compatible types
    if (dynamic_cast<T*>(rhs.GetPointerValue())) {
        Link(rhs);
    }
    return rhs;
}

template <class T>
template <class U>
Cube::Handle<U>& Cube::Handle<T>::operator = (Cube::Handle<U>& rhs) {
    // Going to replace the value of this smart pointer, so unref and
    // possible delete.
    Unlink();
    // Compatible types
    if (dynamic_cast<T*>(rhs.GetPointerValue())) {
        Link(rhs);
    }
    return rhs;
}

template <class T>
template <class U>
const Cube::Handle<U>& Cube::Handle<T>::operator = (const Cube::Handle<U>& rhs) {
    // Going to replace the value of this smart pointer, so unref and
    // possible delete.
    Unlink();
    // Compatible types
    if (dynamic_cast<T*>(rhs.GetPointerValue())) {
        Link(rhs);
    }
    return rhs;
}

template <class T>
T& Cube::Handle<T>::operator*() const {
    TObject* object = GetPointerValue();
    if (!object) {
        CUBE_ERROR << "Dereferencing a NULL handle " << typeid(T).name()
                   << std::endl;
        throw std::runtime_error("Bad reference");
    }
    T* pointer = dynamic_cast<T*>(object);
    if (!pointer) {
        CUBE_ERROR << "Dereferencing with an invalid cast "
                   << typeid(T).name()
                   << std::endl;
        throw std::runtime_error("Bad reference");
    }
    return *pointer;
}

template <class T>
T* Cube::Handle<T>::operator->() const {
    TObject* object = GetPointerValue();
    if (!object) {
        CUBE_ERROR << "Referencing a NULL handle " << typeid(T).name()
                   << std::endl;
        throw std::runtime_error("Bad reference");
    }
    T* pointer = dynamic_cast<T*>(object);
    if (!pointer) {
        CUBE_ERROR << "Referencing with an invalid cast"
                   << typeid(T).name()
                   << std::endl;
        throw std::runtime_error("Bad reference");
    }
    return pointer;
}
#endif

#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// End:

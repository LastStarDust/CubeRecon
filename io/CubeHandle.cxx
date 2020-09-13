#include <iostream>
#include <map>
#include <set>

#include <TROOT.h>
#include <TClass.h>

#include "CubeHandle.hxx"


namespace {
    int gHandleBaseCount = 0;
    int gLastHandleCount = 0;
}

ClassImp(Cube::HandleBase)
Cube::HandleBase::HandleBase() : fCount(0), fHandleCount(0) {
    ++gHandleBaseCount;
}
Cube::HandleBase::~HandleBase() {
    --gHandleBaseCount;
}

ClassImp(Cube::HandleBaseDeletable)
Cube::HandleBaseDeletable::HandleBaseDeletable()
    : fObject(NULL) { }
Cube::HandleBaseDeletable::HandleBaseDeletable(TObject* pointee)
    : fObject(pointee) { }
Cube::HandleBaseDeletable::~HandleBaseDeletable() {
    DeleteObject();
}
void Cube::HandleBaseDeletable::DeleteObject() {
    if (!fObject) return;
    // Actually delete the object.
    if (IsOwner()) delete fObject;
    fObject = NULL;
}

ClassImp(Cube::HandleBaseUndeletable)
Cube::HandleBaseUndeletable::HandleBaseUndeletable() : fObject(NULL) { }
Cube::HandleBaseUndeletable::HandleBaseUndeletable(TObject* pointee)
    : fObject(pointee) { }
Cube::HandleBaseUndeletable::~HandleBaseUndeletable() {
    DeleteObject();
}
void Cube::HandleBaseUndeletable::DeleteObject() {
    fObject = NULL;  // Just set the object pointer to NULL;
}

bool Cube::CleanHandleRegistry(bool) {
    bool result = (gHandleBaseCount==gLastHandleCount);
    if (!result) {
        CUBE_ERROR << "CleanHandleRegistry::"
                   << " Handle Count: " << gHandleBaseCount
                   << " Change: " << gHandleBaseCount - gLastHandleCount
                   << std::endl;;
        gLastHandleCount = gHandleBaseCount;
    }
    return result;
}

ClassImp(Cube::VHandle)
Cube::VHandle::VHandle() {Default(NULL);}
Cube::VHandle::~VHandle() {}

void Cube::VHandle::Default(Cube::HandleBase* handle) {
    fHandle = handle;
    SetBit(kWeakHandle,false);
    if (fHandle) {
        fHandle->CheckHandle();
        fHandle->IncrementReferenceCount();
        fHandle->IncrementHandleCount();
    }
}

void Cube::VHandle::Link(const Cube::VHandle& rhs) {
    // Copy the handle.
    fHandle = rhs.fHandle;
    if (!fHandle) return;
    fHandle->CheckHandle();
    fHandle->IncrementHandleCount();
    if (IsWeak()) return;
    fHandle->IncrementReferenceCount();
}

void Cube::VHandle::Unlink() {
    if (!fHandle) return;
    fHandle->CheckHandle();
    if (!IsWeak()) fHandle->DecrementReferenceCount();
    fHandle->DecrementHandleCount();
    CheckSurvival();
    fHandle = NULL;
}

void Cube::VHandle::MakeWeak() {
    if (IsWeak()) return;
    SetBit(kWeakHandle,true);
    // Decrement the reference count to the object, but leave the handle count
    // unchanged.
    if (!fHandle) return;
    fHandle->CheckHandle();
    fHandle->DecrementReferenceCount();
    CheckSurvival();
}

void Cube::VHandle::MakeLock() {
    if (!IsWeak()) return;
    SetBit(kWeakHandle,false);
    // Increment the reference count to the object, but leave the handle count
    // unchanged, but only if there is a valid handle, and a valid object.
    if (!fHandle) return;
    if (!fHandle->GetObject()) return;
    fHandle->CheckHandle();
    fHandle->IncrementReferenceCount();
}

void Cube::VHandle::CheckSurvival() {
    // The handle doesn't exist, so just return.
    if (!fHandle) return;
    // Check for old handles.
    fHandle->CheckHandle();
    // The handle counter is zero so nothing (no strong, or weak handles) is
    // using this THandleBase and it should be deleted.  This also deletes the
    // object.
    if (fHandle->GetHandleCount() < 1) {
        fHandle->DeleteObject();
        delete fHandle;
        fHandle = NULL;
        return;
    }
    // The reference counter is zero, so no strong handles are referencing the
    // object.  Delete the object, but leave the THandleBase.
    if (fHandle->GetReferenceCount() < 1) {
        fHandle->DeleteObject();
    }
}

TObject* Cube::VHandle::GetPointerValue() const {
    if (!fHandle) return NULL;
    return fHandle->GetObject();
}

void Cube::VHandle::Release(void) {
    if (!fHandle) return;
    fHandle->Release();
}

bool Cube::VHandle::operator == (const Cube::VHandle& rhs) const {
    if (fHandle == rhs.fHandle) return true;
    if (!fHandle) return false;
    if (!rhs.fHandle) return false;
    return (fHandle->GetObject() == rhs.fHandle->GetObject());
}

void Cube::VHandle::ls(Option_t *opt) const {
    TROOT::IndentLevel();
    std::cout << ClassName() << "(" << this << "):: ";
    if (strstr(opt,"size")) {
        TClass* cls = Class();
        if (!cls) return;
        std::cout << " (" << cls->Size() << " b)";
    }
    if (fHandle) {
        std::cout << " Refs: " << fHandle->GetReferenceCount();
        std::cout << " (" << fHandle->GetHandleCount() << ")";
        if (IsWeak()) std::cout << " (weak)";
        else if (fHandle->IsOwner()) std::cout << " (owner)";
        else std::cout << " (released)";
    }
    std::cout << std::endl;
    TROOT::IncreaseDirLevel();
    const TObject* ptr = GetPointerValue();
    if (ptr) ptr->ls(opt);
    TROOT::DecreaseDirLevel();
}

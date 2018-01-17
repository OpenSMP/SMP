//
// Created by Lu WJ on 10/04/2017
//

#ifndef CRYPTCONV_COWPTR_HPP
#define CRYPTCONV_COWPTR_HPP
#include <memory>

template <class T>
class CowPtr
{
public:
    typedef std::shared_ptr<T> RefPtr;

private:
    RefPtr m_sp;

    void detach() {
        T* tmp = m_sp.get();
        if( !( tmp == 0 || m_sp.unique() ) ) {
            m_sp = RefPtr( new T( *tmp ) );
        }
    }

public:
    CowPtr(T* t = nullptr) : m_sp(t) {}

    CowPtr(const RefPtr& refptr) : m_sp(refptr) {}

    const T& operator*() const { return *m_sp; }

    T& operator*() { detach(); return *m_sp; }

    const T* operator->() const { return m_sp.operator->(); }

    T* operator->() { detach(); return m_sp.operator->(); }
};

#endif //CRYPTCONV_COWPTR_HPP

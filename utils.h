/*  utils.h - Assorted template factory, singleton factory and dispatcher class declarations.

     Mike Jackson, The University of Edinburgh & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

    Copyright (C) 2015 University of Oxford  */


#ifndef __FABBER_UTILS_H
#define __FABBER_UTILS_H 1

#include <map>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

/**
 * Template factory class.
 *
 * Manages a map from names to function pointers. Each function
 * pointer is assumed to point to a function that returns a pointer to
 * an object of type T.
 */
template <class T>
class TemplateFactory {
  public:
    /**
     * Function pointer to a zero argument function that returns a
     * pointer to an object of type T.
     */
    typedef T* (*Function)(void);
    /** Constructor. */
    TemplateFactory();
    /** Destructor. */
    ~TemplateFactory();
    /**
     * Add a mapping from a given name to a function pointer.
     * @param name Name by which the function pointer can be
     * accessed.
     * @param function Function pointer.
     */
    void Add(const string& name, Function function);
    /**
     * Invoke the function pointer with the given name and return a
     * pointer to an object.
     * @param name
     * @return pointer, or NULL if name is
     * not known.
     */ 
    T* Create(const string& name);
    /** 
     * Get the list of names.
     * @return list of names.
     */
    vector<string> GetNames();
    /** 
     * Is there a function pointer with the given name?
     * @param name
     * @return true if so, false otherwise.
     */
    bool HasName(const string& name);
  private:
    /** Map from names to function pointers. */
    map<string, Function> functionMap_;
};

template <class T>
TemplateFactory<T>::TemplateFactory() {
  // No-op.
}

template <class T>
TemplateFactory<T>::~TemplateFactory() {
  functionMap_.clear(); 
}

template <class T>
vector<string> TemplateFactory<T>::GetNames() {
  vector<string> names;
  for (typename map<string, Function>::iterator it = 
       functionMap_.begin(); it != functionMap_.end(); it++) {
    names.push_back(it->first);
  }
  return names;
}

template <class T>
bool TemplateFactory<T>::HasName(const string& name) {
  return (functionMap_.find(name) != functionMap_.end());
}

template <class T>
void TemplateFactory<T>::Add(const string& name, Function function) {
  functionMap_[name] = function;
}

template <class T>
T* TemplateFactory<T>::Create(const string& name) {
  if (functionMap_.count(name)) {
    return functionMap_[name]();
  }
  return NULL;
}



/**
 * Singleton template factory class.
 *
 * Maintains a singleton instance of a \ref TemplateFactory.
 */
template <class T>
class SingletonFactory : public TemplateFactory<T> {
  public:
    /**
     * Returns pointer to singleton instance of this class.
     * @return instance.
     */
    static SingletonFactory* GetInstance();
    /**
     * Delete the singleton instance.
     */
    static void Destroy();
  private:
    /** Singleton instance of this class. */
    static SingletonFactory* singleton_;
    /** Constructor. */
    SingletonFactory();
};

template <class T>
SingletonFactory<T>::SingletonFactory() {
  // No-op.
}

template <class T>
void SingletonFactory<T>::Destroy() {
  if (singleton_ != NULL) {
    delete singleton_;
    singleton_ = NULL;
  }
}

template <class T>
SingletonFactory<T>* SingletonFactory<T>::singleton_ = NULL;

template <class T>
SingletonFactory<T>* SingletonFactory<T>::GetInstance() {
  if (singleton_ == NULL) {
    singleton_ = new SingletonFactory<T>();
  }
  return singleton_;
}



/**
 * Template class for registration of classes with factories.
 * Assumes T supports a GetInstance function which returns an object
 * which supports an Add(name, Function) function where Function
 * is a function pointer.
 * Assumes U supports a NewInstance function that conforms to the type
 * suppoorted by T's Add function.
 */
template<class T, class U>
class FactoryRegistration {
  public:
    /**
     * Constructor.
     * @param name Name by which U's function will be registered,
     */
    FactoryRegistration(string name) {
      (T::GetInstance())->Add(name, &U::NewInstance);
    }
};



/**
 * Simple dispatcher class.
 * 
 * Manages a map from names to function pointers. Each function
 * pointer is assumed to point to a function that carries out some
 * action.
 */
class Dispatcher {
  public:
    /**
     * Function pointer to a zero argument function that does some
     * action, but returns no result.
     */
    typedef void (*Function)(void);
    /** Constructor. */
    Dispatcher();
    /** Destructor. */
    ~Dispatcher();
    /**
     * Add a mapping from a given name to a function pointer.
     * @param name Name by which the function pointer can be
     * accessed.
     * @param function Function pointer.
     */
    void Add(const string& name, Function function);
    /**
     * Invoke the function pointer with the given name.
     * @param name
     */ 
    void Dispatch(const string& name);
    /** 
     * Get the list of names.
     * @return list of names.
     */
    vector<string> GetNames();
    /** 
     * Is there a function pointer with the given name?
     * @param name
     * @return true if so, false otherwise.
     */
    bool HasName(const string& name);
  private:
    /** Map from names to function pointers. */
    map<string, Function> functionMap_;
};

#endif // __FABBER_UTILS_H

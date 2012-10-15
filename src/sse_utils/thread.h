/**
    libgapmis: a library for pairwise sequence aligment with a single gap.
    Copyright (C) 2012 Nikolaos Alachiotis, Simon Berger, Tomas Flouri, and
    Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

// this is a linux only copy of the thread library from github.com/sim82/ivy_mike.
// it is a 'drop-in' replacement for either boost::thread or std::thread from c++11


#ifndef __ivy_mike__thread_h
#define __ivy_mike__thread_h


#include <pthread.h>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cerrno>



class thread {
    
public:
    typedef pthread_t native_handle_type;
private:
    native_handle_type m_thread;
    bool m_valid_thread;
    
    thread( const thread &other );
    const thread &operator=( const thread &other );
    
    template<typename callable>
    static void *call( void *f ) {
        std::auto_ptr<callable> c(static_cast<callable *>(f));
        
         try {
            (*c)();
        } catch( std::runtime_error x ) {
            std::cerr << "uncaught std::runtime_error in ivy_mike::thread: " << x.what() << std::endl;
//             std::cerr << x.what() << std::endl;
            
            throw;
        } catch( std::exception x ) {
            std::cerr << "uncaught std::exception in ivy_mike::thread" << std::endl;
//             std::cerr << x.what() << std::endl;
            
            throw;
        }
       
        return 0;
    }
    
public:
    
    thread() : m_valid_thread( false ) {
        
    }
    
    template<typename callable>
    thread( const callable &c ) {
        // dynamically allocating a copy of callable seems to be the only way to get the object into the 
        // thread start function without specializing the whole thread object...
        int ret = pthread_create( &m_thread, 0, call<callable>, new callable(c) );
        
        if( ret == EAGAIN ) {
            throw std::runtime_error( "could not create thread: resource_unavailable_try_again" );
        } else if( ret != 0 ) {
            throw std::runtime_error( "could not create thread: internal error" );
        } else {
            m_valid_thread = true;
        }
    }
    
    ~thread() {
        // this (=nothing) is actually the right thing to do (TM) if I read 30.3.1.3 if the iso c++0x standard corrctly
        // ^^^ WTF? shouldn't this cause an ugly abort?
        
//         std::cout << "~thread\n"; 
        if( joinable() ) {
            std::cerr << "ivy_mike::thread warning: destructor of joinable thread called. possible ressource leak" << std::endl;
        }
    }
    
    void swap( thread &other ) {
        std::swap( m_thread, other.m_thread );
        std::swap( m_valid_thread, other.m_valid_thread );
    }
    
    bool joinable() {
        return m_valid_thread;
    }
    
    void join() {
        
        if( joinable() ) {
            void *rv;
            pthread_join(m_thread, &rv );
            m_valid_thread = false;
        }
    }
    
    native_handle_type native_handle() {
        if( !joinable() ) {
            throw std::runtime_error( "native_handle: thread not joinable" );
        }
        return m_thread;
    }
};



class mutex {
    pthread_mutex_t m_mtx;

public:
    
    mutex() {
        pthread_mutex_init( &m_mtx, 0);
    }
    ~mutex() {
        pthread_mutex_destroy(&m_mtx);
    }
    
    inline void lock() {
        pthread_mutex_lock(&m_mtx);
    }
    
    inline void unlock() {
        pthread_mutex_unlock(&m_mtx);
    }
};



class thread_group {
    std::vector<thread *> m_threads;
    
public:
    
    ~thread_group() {
//         std::cout << "thread_group destructor: fallback join:\n";
//         for( std::vector<thread *>::iterator it = m_threads.begin(); it != m_threads.end(); ++it ) {
//             std::cout << "joinable: " << (*it)->joinable() << "\n";
//         }
        
        join_all();
        
        for( std::vector<thread *>::iterator it = m_threads.begin(); it != m_threads.end(); ++it ) {
            delete *it;
        }
        
    }
    
    
    
    template<typename callable>
    void create_thread( const callable &c ) {
        m_threads.push_back(0); // may throw. so pre-allocate before the thread is created
        m_threads.back() = new thread( c );
        
        
    }
    
    
    void join_all() throw() {
        
        try {
            for( std::vector<thread *>::iterator it = m_threads.begin(); it != m_threads.end(); ++it ) {
                if( (*it)->joinable() ) {
                    (*it)->join();
                }
                
                //delete (*it);
            }
            
            m_threads.resize(0);
        } catch(...) {
            std::cerr << "BUG: unexpected exception in thread_group::join_all\n"; // kind of stupid: printing this message might throw...
        }
    }
    
    size_t size() {
        return m_threads.size();
    }
};

inline void swap( thread &t1, thread &t2 ) {
    t1.swap(t2);
}


template <typename mtx_t>
class lock_guard {
    mtx_t &m_mtx;
    

    lock_guard( const lock_guard & );
    lock_guard &operator=(const lock_guard & );

public:
    lock_guard( mtx_t &mtx ) 
     : m_mtx(mtx)
    {
        m_mtx.lock();
    }
    
    ~lock_guard() {
        m_mtx.unlock();
    }
};

#endif

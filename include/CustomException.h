/*
 * MyException.hpp
 *
 *  Created on: 03.09.2016
 *      Author: Ruhullah
 */

#ifndef CUSTOMEXCEPTION_H_
#define CUSTOMEXCEPTION_H_

#include <exception>
#include <string>

//! The Exception class from which every exception type in this project inherits from.
/*!
 * This is the superclass for all the custom exception types used in this project.
 * This class itself inherits from std::exception.
 */
class CustomException: public std::exception {
public:
    /*!
     * The constructor taking a message as parameter.
     * @param msg The message of the exception.
     */
    CustomException(const std::string& msg) :
            msg(msg)
    {
    }

    /*!
     * The virtual destructor.
     */
    virtual ~CustomException() throw ()
    {
    }

    /*!
     * Returns the message of the exception.
     * @return
     */
    virtual const char* what() const throw ()
    {
        return msg.c_str();
    }

private:
    std::string msg;
};

#endif /* CUSTOMEXCEPTION_H_ */

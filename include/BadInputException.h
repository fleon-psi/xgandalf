/*
 * BadInput.hpp
 *
 *  Created on: 03.09.2016
 *      Author: Ruhullah
 */

#ifndef BADINPUTEXCEPTION_H_
#define BADINPUTEXCEPTION_H_

#include <CustomException.h>

//! A specialization of MyException.
/*!
 * This exception is thrown whenever the input to a method/function is not as expected
 */
class BadInputException: public CustomException {
public:
    BadInputException(const std::string& msg) :
            CustomException(msg)
    {
    }

    virtual ~BadInputException() throw ()
    {
    }
};

#endif /* BADINPUTEXCEPTION_H_ */

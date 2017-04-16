/*
 * WrongUsageException.h
 *
 *  Created on: 16.04.2017
 *      Author: Yaro
 */

#ifndef WRONGUSAGEEXCEPTION_H_
#define WRONGUSAGEEXCEPTION_H_

#include <CustomException.h>

//! A specialization of MyException.
/*!
 * This exception is thrown whenever a function/method is used in a wrong way
 */
class WrongUsageException: public CustomException {
public:
    WrongUsageException(const std::string& msg) :
            CustomException(msg)
    {
    }

    virtual ~WrongUsageException() throw ()
    {
    }
};

#endif /* WRONGUSAGEEXCEPTION_H_ */

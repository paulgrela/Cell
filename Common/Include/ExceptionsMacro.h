#pragma once

#ifndef _EXCEPTIONS_MACRO_H_
#define _EXCEPTIONS_MACRO_H_

#include <iostream>
#include <exception>

#include "Logger.h"


#define CATCH(str) \
		catch (const std::exception& ExceptionObject) \
		{ \
			LoggersManagerObject.Log(STREAM("Exception caught in " << str << ": " << ExceptionObject.what())); \
		} \
		catch (...) \
		{ \
			LoggersManagerObject.Log(STREAM("Unknown exception caught in " << str));  \
		}

#define CATCH_AND_THROW(str) \
		catch (const std::exception& ExceptionObject) \
		{ \
			LoggersManagerObject.Log(STREAM("Exception caught in " << str << ": " << ExceptionObject.what())); \
			throw; \
		} \
		catch (...) \
		{ \
			LoggersManagerObject.Log(STREAM("Unknown exception caught in " << str));  \
			throw; \
		} 

#define CATCH_AND_WORK(str, Before, After) \
		catch (const std::exception& ExceptionObject) \
		{ \
		    Before; \
			LoggersManagerObject.Log(STREAM("Exception caught in " << str << ": " << ExceptionObject.what())); \
			After; \
		} \
		catch (...) \
		{ \
		    Before; \
			LoggersManagerObject.Log(STREAM("Unknown exception caught in " << str));  \
			After; \
		}


#define CATCH_COUT(str) \
		catch (const std::exception& ExceptionObject) \
		{ \
			std::cout << "Exception caught in " << str << ": " << ExceptionObject.what() << endl; \
		} \
		catch (...) \
		{ \
			std::cout << "Unknown exception caught in " << str << endl;  \
		} 

#define CATCH_AND_THROW_COUT(str) \
		catch (const std::exception& ExceptionObject) \
		{ \
			std::cout << "Exception caught in " << str << ": " << ExceptionObject.what() << endl; \
			throw; \
		} \
		catch (...) \
		{ \
			std::cout << "Unknown exception caught in " << str << endl;  \
			throw; \
		} 

#define CATCH_AND_WORK_COUT(str, Before, After) \
		catch (const std::exception& ExceptionObject) \
		{ \
		    Before; \
			std::cout << "Exception caught in " << str << ": " << ExceptionObject.what() << endl; \
			After; \
		} \
		catch (...) \
		{ \
		    Before; \
			std::cout << "Unknown exception caught in " << str << endl;  \
			After; \
		}

#endif
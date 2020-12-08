
#ifndef Log_hpp
#define Log_hpp

#ifdef USE_BOOST_LOG
	#include <boost/log/core.hpp>
	#include <boost/log/trivial.hpp>
	#include <boost/log/expressions.hpp>
	#include <boost/log/sinks/text_file_backend.hpp>
	#include <boost/log/utility/setup/file.hpp>
	#include <boost/log/utility/setup/common_attributes.hpp>
	#include <boost/log/sources/severity_logger.hpp>
	#include <boost/log/sources/record_ostream.hpp>
	#include <boost/log/utility/setup/settings.hpp>
	#include <boost/log/utility/setup/from_settings.hpp>
	#include <boost/date_time/posix_time/posix_time.hpp>

	#include <boost/format.hpp>

	inline void initLog(int rank,int severity){

		std::string file=str(boost::format("log.%d") % rank);
		boost::log::settings setts;

		setts["Core"]["Disablelog"] = false;
		setts["Sinks.Destination"] = "TextFile";
		setts["Sinks.File.Destination"] = "TextFile";
		setts["Sinks.File.FileName"] = file;
		setts["Sinks.File.AutoFlush"] = true;
		setts["Sinks.File.Format"] = "[ %TimeStamp% ] %Severity%: %Message%";

		//keywords::format = "[%TimeStamp%]: %Message%"
		boost::log::register_simple_formatter_factory< boost::log::trivial::severity_level, char >("Severity");

		boost::log::add_common_attributes();

		boost::log::init_from_settings(setts);

		switch(severity) {
		case 6:
			boost::log::core::get()->set_filter
			(
				boost::log::trivial::severity >= boost::log::trivial::fatal
			);
			break;
		case 5:
			boost::log::core::get()->set_filter
			(
				boost::log::trivial::severity >= boost::log::trivial::error
			);
			break;
		case 4:
			boost::log::core::get()->set_filter
			(
				boost::log::trivial::severity >= boost::log::trivial::warning
			);
			break;
		case 3:
			boost::log::core::get()->set_filter
			(
				boost::log::trivial::severity >= boost::log::trivial::info
			);
			break;
		case 2:
			boost::log::core::get()->set_filter
			(
				boost::log::trivial::severity >= boost::log::trivial::debug
			);
			break;
		case 1:
			boost::log::core::get()->set_filter
			(
				boost::log::trivial::severity >= boost::log::trivial::trace
			);
			break;
		default:
			//setts["Core"]["Disablelog"] = true;
			break;
		}


		boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;


	};

	#ifdef VERBOSE
		#define LOGGER(str) BOOST_LOG_SEV(BaseEngine::lg, boost::log::trivial::trace)<<str;
	#else
		#define LOGGER(str) ;
	#endif

	#define LOGGERA(str) BOOST_LOG_SEV(BaseEngine::lg, boost::log::trivial::trace)<<str;

#else

	#ifdef VERBOSE
		#define LOGGER(str) std::cout<<str<<"\n";
	#else
		#define LOGGER(str) ;
	#endif

	#define LOGGERA(str) std::cout<<str<<"\n";

#endif

#endif

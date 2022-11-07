#include "geogram/basic/assert.h"

namespace GEO
{
	void geo_assertion_failed( const std::string& condition_string, const std::string& file, int line )
	{
		// __debugbreak();
	}
	void geo_should_not_have_reached( const std::string& file, int line )
	{
		__debugbreak();
	}
	AssertMode assert_mode()
	{
		return AssertMode::ASSERT_THROW;
	}
	void set_assert_mode( AssertMode mode )
	{
	}
}  // namespace GEO

#include "geogram/basic/logger.h"

namespace GEO
{
	static std::ostream& getEmptyStream()
	{
		class NullBuffer : public std::streambuf
		{
		  public:
			int overflow( int c )
			{
				return c;
			}
		};
		static NullBuffer null_buffer;
		static std::ostream null_stream( &null_buffer );
		return null_stream;
	}

	std::ostream& Logger::err( const std::string& feature )
	{
		return getEmptyStream();
	}
	std::ostream& Logger::warn( const std::string& feature )
	{
		return getEmptyStream();
	}
	std::ostream& Logger::out( const std::string& feature )
	{
		return getEmptyStream();
	}
	bool Logger::is_initialized()
	{
		return true;
	}
}  // namespace GEO

#include "geogram/basic/command_line.h"

namespace GEO
{
	std::string CmdLine::get_arg( const std::string& name )
	{
		return "";
	}

	bool CmdLine::get_arg_bool( const std::string& name )
	{
		return false;
	}
}  // namespace GEO

#include "geogram/bibliography/bibliography.h"

namespace GEO
{
	namespace Biblio
	{
		void cite( const char* ref, const char* file, int line, const char* function, const char* info )
		{
		}
	}  // namespace Biblio
}  // namespace GEO

#include "geogram/basic/factory.h"

namespace GEO
{
	void InstanceRepo::add( const std::string& name, Instance* instance )
	{
	}
	InstanceRepo::Instance* InstanceRepo::get( const std::string& name )
	{
		return nullptr;
	}
}  // namespace GEO

#include "geogram/basic/numeric.h"

namespace GEO
{
	namespace Numeric
	{
		int32 random_int32()
		{
			return rand();
		}
	}  // namespace Numeric
}  // namespace GEO
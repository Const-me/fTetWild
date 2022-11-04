#include <geogram/basic/environment.h>

namespace GEO
{
	Environment* Environment::instance()
	{
		return nullptr;
	}
	bool Environment::add_environment( Environment* env )
	{
		return false;
	}
	bool Environment::has_value( const std::string& name ) const
	{
		return false;
	}
	bool Environment::get_value( const std::string& name, std::string& value ) const
	{
		return false;
	}
	bool Environment::set_value( const std::string& name, const std::string& value )
	{
		return false;
	}
	bool Environment::add_observer( const std::string& name, VariableObserver* observer )
	{
		return false;
	}
	bool Environment::remove_observer( const std::string& name, VariableObserver* observer )
	{
		return false;
	}
	Environment* Environment::find_environment( const std::string& name )
	{
		return nullptr;
	}
	Environment::~Environment()
	{
	}
	bool Environment::notify_observers( const std::string& name, bool recursive )
	{
		return false;
	}
}  // namespace GEO
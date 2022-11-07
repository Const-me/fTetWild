#pragma once

// Utility class to measure time, and print to console time in seconds
class Timer
{
	const uint64_t started;
	const char* const message;

  public:
	Timer( const char* what );
	~Timer();
	Timer( const Timer& ) = delete;
};
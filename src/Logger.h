/* Logging utility for cossolve.
 *
 * Copyright (C) 2023 Kyle Manke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef COSSOLVE_LOGGER_H
#define COSSOLVE_LOGGER_H

#include <iostream>
#include <chrono>
#include <string_view>
#include <sstream>
#include <iomanip>

namespace cossolve {

enum class LogLevel : uint8_t
{
    debug = 0,
    info = 1,
    warning = 2,
    error = 3,
    fatal = 4
};

template <LogLevel minLevel = LogLevel::info>
class Logger
{
    template <LogLevel level>
    class Entry;

public:

    Logger(std::ostream& out) : out(out) { }

    // Begins a timer for performance measurement
    void startTimer(std::string&& label)
    {
	timerLabel = std::move(label);
	startTime = std::chrono::high_resolution_clock::now();
    }
    // Ends the timer started by startTimer
    void endTimer()
    {
	endTime = std::chrono::high_resolution_clock::now();
	auto us = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
	log() << timerPrefix << timerLabel << timerPostfix
	      << us.count() << timerUnits;
    }

    template <LogLevel level = LogLevel::debug>
    Entry<level> log()
    {
	return Entry<level>(out);
    }

private:
    std::ostream& out;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> endTime;
    std::string timerLabel;

    static constexpr std::string_view logLevelLabels[] =
    {
	" [DEBUG] ",
	" [INFO] ",
	" [WARNING] ",
	" [ERROR] ",
	" [FATAL] "
    };
    static constexpr std::string_view timerPrefix = "Operation `";
    static constexpr std::string_view timerPostfix = "` completed in ";
    static constexpr std::string_view timerUnits = " us.";

    template <LogLevel level>
    class Entry
    {
    public:
	// Constructor generates the timestamp and stores a reference to the output stream
        Entry(std::ostream& out) : out(out)
	{
	    if constexpr(level >= minLevel)
	    {
		std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		auto timeStamp = std::put_time(std::localtime(&now), "[%Y-%m-%d %X]");
		str << timeStamp << Logger::logLevelLabels[static_cast<uint8_t>(level)];
	    }
	    return;
	}
	// Upon destruction, output `str` to the output stream
	~Entry()
	{
	    if constexpr(level >= minLevel)
	    {
		str << std::endl;
		out << str.str();
	    }

	    return;
	}
	// Forward << to the underlying stream
	template <typename T>
	decltype(auto) operator<< (T&& t)
	{
	    if constexpr(level >= minLevel)
	    {
		return str << std::forward<T>(t);
	    }
	    else
	    {
		return *this;
	    }
	}
	
    private:
	std::stringstream str;
	std::ostream& out;
    };
};
    
} // namespace cossolve

#endif // COSSOLVE_LOGGER_H

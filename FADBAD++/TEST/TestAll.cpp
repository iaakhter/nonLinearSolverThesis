#include "STDReportLog.h"

#include "TestFADBAD.h"
#include "TestTAD.h"

int main()
{
	STDReportLog log(std::cerr);

	TestFADBAD testFADBAD;
	testFADBAD.run(log);

	TestTAD testTAD;
	testTAD.run(log);

	return 0;
}

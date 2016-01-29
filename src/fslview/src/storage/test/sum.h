// sum.h: interface for the Sum class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SUM_H)
#define SUM_H

#include "../volumevisitor.h"

class Sum: public VolumeVisitor
{
public:
	Sum();

	virtual void visit(VolumeB::WeakHandle target)
	{
		calculate(target);
	}
	virtual void visit(VolumeUS::WeakHandle target)
	{
		calculate(target);
	}
	virtual void visit(VolumeF::WeakHandle target)
	{
		calculate(target);
	}

	float total() const { return total_; }
	unsigned int count() const { return count_; }

private:

	template <typename VoxelType> 
		void calculate(VolumeStore<VoxelType>::WeakHandle target);
				
private:
	float total_;
	unsigned int count_;
};

#endif 

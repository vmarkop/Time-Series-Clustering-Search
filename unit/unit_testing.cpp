#include "mathUtils.h"
#include "projectUtils.h"

#include "gtest/gtest.h"


/*  Test euclideanDistance()
 *  when a NULL point is given.
 */
TEST(euclideanDistance, null_point)
{
    PointPtr point = new Point__Struct;
    point->id = "id";
    point->coords.push_back(1.1);
    point->coords.push_back(2.2);

    ASSERT_EQ(-1,euclideanDistance(point,NULL ,2));
    ASSERT_EQ(-1,euclideanDistance(NULL ,point,2));
    ASSERT_EQ(-1,euclideanDistance(NULL ,NULL,20));

    delete point;
}

/*  Test euclideanDistance()
 *  when the same point is
 *  given in both arguments.
 */
TEST(euclideanDistance, same_point)
{
    PointPtr point = new Point__Struct;
    point->id = "id";
    point->coords.push_back(1.1);
    point->coords.push_back(2.2);

    ASSERT_EQ(0,euclideanDistance(point,point,2));

    delete point;
}

/*  Test giving buildTree()
 *  function 0 height.
 */
TEST(buildTree, zeroheight)
{
    ASSERT_EQ(NULL,buildTree(0));
}

/*  Showcase correct
 *  behaviour of
 *  function ^2
 */
TEST(powerWithBase2, usage)
{
    ASSERT_EQ(8,powerWithBase2(3));
    ASSERT_EQ(1,powerWithBase2(0));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
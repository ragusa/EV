ADD_EXECUTABLE(group_fe_values group_fe_values.cc)
DEAL_II_SETUP_TARGET(group_fe_values)
SET(TEST_TARGET group_fe_values)
DEAL_II_PICKUP_TESTS()

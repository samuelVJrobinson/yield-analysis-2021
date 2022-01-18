# Yield simulations given results from model:
# i.e. "What happens to overall yield if I change boundary X to Y? (Given field size of ...)"

#Requires Galpern lab machine

#Idea: 
# - Get cells in the shape of a quarter section 
# - Calculate boundary distance
# - Simulate yield N times, get average for field, save
# - Repeat for different boundary types:
#   - Standard-Shelterbelt combinations
#   - Standard-Wetland combination
#   - Shelterbelt-Wetland combinations
# - Repeat for canola, peas, wheat
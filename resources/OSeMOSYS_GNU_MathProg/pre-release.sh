string=$1
export prerelease='false'
if [[ $string == *"alpha"* ]]; then
  prerelease='true'
elif [[ $string == *"beta"* ]]; then
  prerelease='true'
else
  prerelease='false'
fi
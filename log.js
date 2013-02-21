/*
  Function that prints log messages.
  

  Copyright 2013 Guillaume Lathoud
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
  A copy of the Apache License Version 2.0 as of February 20th, 2013
  can be found in the file ./LICENSE.TXT
*/

(function (G) {

    if (G.log)
        return;

    G.log = log;

    function log(s)
    {
        if (typeof console !== 'undefined')
            console.log(s);  // Browsers
        else
            print(s);  // V8, Rhino
    }

})(this);

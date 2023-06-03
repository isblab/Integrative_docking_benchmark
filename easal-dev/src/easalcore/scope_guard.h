/*
 This file is part of EASAL.

 EASAL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 EASAL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SCOPE_GUARD_H
#define SCOPE_GUARD_H
template <class Fun>
class scope_guard {
 public:
  scope_guard(Fun f) : fun_(std::move(f)), enabled_(true) {}

  scope_guard(scope_guard&& x) : fun_(std::move(x.fun_)), enabled_(x.enabled_) {
    x.enabled_ = false;
  }

  ~scope_guard() {
    if (enabled_) fun_();
  }

 private:
  Fun fun_;
  bool enabled_;
};

// Creates a guard that executes `f` as soon as it goes out of scope.
template <class Fun>
scope_guard<Fun> make_scope_guard(Fun f) {
  return {std::move(f)};
}

#endif

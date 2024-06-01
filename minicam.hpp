/********************************************************************************
BSD 3-Clause License

Copyright (c) 2023, Li Yunqiang (walkfish8@hotmail.com)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/
#pragma once

#ifndef _MINICAM_MINICAM_HPP_
#define _MINICAM_MINICAM_HPP_

#include <string>
#include <vector>
#include <iostream>

namespace minicam
{

template <typename T> struct Point2_ {
    T x, y;
};
template <typename T> struct Point3_ {
    T x, y, z;
};

/// @brief 基础数据类型为数值类型
template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
struct CoordType {
    typedef T type;
};
/// @brief 辅助函数，查询类的基础数据类型;
template <typename T> struct CoordTraits: CoordType<T> {};
template <typename T> using coord_traits_t = typename CoordTraits<T>::type;

enum class EulerAxisMode : int { XYZ = 0, XZY, YXZ, YZX, ZXY, ZYX };

/// @brief SO3旋转类
/// @tparam Float 浮点类型
template <typename Float> struct Rotation_ {
    static_assert(std::is_floating_point<Float>::value,
        "Struct Rotation_ Can Only Be Instantiated With Float Point Type.");
    /// @brief 初始化，为单位矩阵.
    Rotation_() : a1(1), a2(0), a3(0), b1(0), b2(1), b3(0), c1(0), c2(0), c3(1)
    {}
    /// @brief 通过3x3的矩阵进行初始化;
    Rotation_(Float data[9]) { memcpy(&a1, data, sizeof(Rotation_)); }
    Rotation_(const Rotation_& R) { memcpy(&a1, &R, sizeof(Rotation_)); }
    /// @brief 通过欧拉角初始化,
    /// @param euler_rad 旋转弧度
    /// @param axis [0,1,2]分别对应axis[x,y,z]三个轴.
    Rotation_(Float euler_rad, int axis = 0);
    /// @brief 通过旋转向量初始化
    Rotation_(Float rx, Float ry, Float rz);
    /// @brief 通过四元素初始化
    Rotation_(Float qw, Float qx, Float qy, Float qz);
    Rotation_(Float m1, Float m2, Float m3, Float m4, Float m5, Float m6,
        Float m7, Float m8, Float m9)
        : a1(m1), a2(m2), a3(m3), b1(m4), b2(m5), b3(m6), c1(m7), c2(m8), c3(m9)
    {}
    // RULERMVS_SMART_CONVERT_MEMBER_FUNC(Rotation_<Float>)

    Rotation_ t() const { return {a1, b1, c1, a2, b2, c2, a3, b3, c3}; }
    Rotation_ inv() const
    {
        // 这里实现了三维矩阵的逆，不进行判断；
        Float det = a1 * (b2 * c3 - b3 * c2) + b1 * (c2 * a3 - a2 * c3) +
                    c1 * (a2 * b3 - b2 * a3);
        Float A11 = (b2 * c3 - b3 * c2) / det;
        Float A12 = (b3 * c1 - b1 * c3) / det;
        Float A13 = (b1 * c2 - b2 * c1) / det;
        Float A21 = (a3 * c2 - a2 * c3) / det;
        Float A22 = (a1 * c3 - a3 * c1) / det;
        Float A23 = (a2 * c1 - a1 * c2) / det;
        Float A31 = (a2 * b3 - a3 * b2) / det;
        Float A32 = (a3 * b1 - a1 * b3) / det;
        Float A33 = (a1 * b2 - a2 * b1) / det;
        return {A11, A21, A31, A12, A22, A32, A13, A23, A33};
    }
    Rotation_ dot(const Rotation_& R) const
    {
        Float data[9] = {a1 * R.a1 + a2 * R.b1 + a3 * R.c1,
            a1 * R.a2 + a2 * R.b2 + a3 * R.c2,
            a1 * R.a3 + a2 * R.b3 + a3 * R.c3,
            b1 * R.a1 + b2 * R.b1 + b3 * R.c1,
            b1 * R.a2 + b2 * R.b2 + b3 * R.c2,
            b1 * R.a3 + b2 * R.b3 + b3 * R.c3,
            c1 * R.a1 + c2 * R.b1 + c3 * R.c1,
            c1 * R.a2 + c2 * R.b2 + c3 * R.c2,
            c1 * R.a3 + c2 * R.b3 + c3 * R.c3};
        return data;
    }
    template <typename Tp> Point3_<Tp> rotate(const Point3_<Tp>& pt) const
    {
        return {(Tp)(a1 * pt.x + a2 * pt.y + a3 * pt.z),
            (Tp)(b1 * pt.x + b2 * pt.y + b3 * pt.z),
            (Tp)(c1 * pt.x + c2 * pt.y + c3 * pt.z)};
    }
    template <typename Tp>
    void rotate(const Point3_<Tp>* src, Point3_<Tp>* dst, size_t pt_num) const
    {
        MVS_OMP_PARALLEL_FOR
        for (int i = 0; i < (int)pt_num; ++i) dst[i] = rotate(src[i]);
    }
#ifdef RULERMVS_USE_SSE
    void rotate(const Point3f* src, Point3f* dst, size_t num) const;
#endif
    template <typename Tp> void rotate(const std::vector<Point3_<Tp>>& src,
        std::vector<Point3_<Tp>>& dst) const
    {
        if (src.empty()) return;
        if (src.size() != dst.size()) dst.resize(src.size());
        rotate<Tp>(&src[0], &dst[0], src.size());
    }
    Rotation_& operator*=(const Rotation_& R) { return *this = dot(R); }

    friend inline std::ostream& operator<<(std::ostream& O, const Rotation_& R)
    {
        return O << "Rotation(" << R.a1 << "," << R.a2 << "," << R.a3 << ","
                 << R.b1 << "," << R.b2 << "," << R.b3 << "," << R.c1 << ","
                 << R.c2 << "," << R.c3 << ")";
    };

    Float a1, a2, a3, b1, b2, b3, c1, c2, c3;
};
using Rotation = Rotation_<double>;
template <typename T> struct CoordTraits<Rotation_<T>>: CoordType<T> {};
template <typename Float>
static inline Rotation_<Float> operator-(const Rotation_<Float>& R)
{
    return {-R.a1, -R.a2, -R.a3, -R.b1, -R.b2, -R.b3, -R.c1, -R.c2, -R.c3};
}
template <typename Float> static inline Rotation_<Float> operator*(
    const Rotation_<Float>& R1, const Rotation_<Float>& R2)
{
    return R1.dot(R2);
}
template <typename Float> static inline bool operator==(
    const Rotation_<Float> R1, const Rotation_<Float>& R2)
{
    const Float e = std::numeric_limits<Float>::epsilon();
    for (int i = 0; i < 9; ++i)
        if (std::abs((&R1.a1)[i] - (&R2.a1)[i]) > e) return false;
    return true;
}

/// @brief  SE3姿态类,用于表达姿态和位置信息.
/// @tparam Float 浮点类型，一般的我们推荐用双精度double类型.
template <typename Float> struct Pose_: Rotation_<Float>, Point3_<Float> {
    static_assert(std::is_floating_point<Float>::value,
        "struct Pose can only be instantiated with float point type.");
    typedef Point3_<Float>   OBase;
    typedef Rotation_<Float> RBase;
    Pose_() : RBase(), OBase(0, 0, 0) {}
    Pose_(const std::string& path) { load(path); }
    Pose_(const Pose_& pose) : RBase(pose), OBase(pose) {}
    Pose_(const RBase& R, const OBase& O) : RBase(R), OBase(O) {}
    Pose_(Float rx, Float ry, Float rz, Float tx, Float ty, Float tz)
        : RBase(rx, ry, rz)
        , OBase(-getRotation().inv().rotate(OBase {tx, ty, tz}))
    {}
    Pose_(Float qw, Float qx, Float qy, Float qz, Float tx, Float ty, Float tz)
        : RBase(qw, qx, qy, qz)
        , OBase(-getRotation().inv().rotate(OBase {tx, ty, tz}))
    {}
    Pose_(Float a1, Float a2, Float a3, Float b1, Float b2, Float b3, Float c1,
        Float c2, Float c3, Float tx, Float ty, Float tz)
        : RBase(a1, a2, a3, b1, b2, b3, c1, c2, c3)
        , OBase(-getRotation().inv().rotate(OBase {tx, ty, tz}))
    {}
    Pose_(const Float mat3x4[12])
        : Pose_(mat3x4[0], mat3x4[1], mat3x4[2], mat3x4[4], mat3x4[5],
              mat3x4[6], mat3x4[8], mat3x4[9], mat3x4[10], mat3x4[3], mat3x4[7],
              mat3x4[11])
    {}
    // RULERMVS_SMART_CONVERT_MEMBER_FUNC(Pose_<Float>)

    bool  load(const std::string& path);
    bool  save(const std::string& path) const;
    OBase center() const { return static_cast<OBase>(*this); }
    RBase getRotation() const { return static_cast<RBase>(*this); }
    OBase getTranslation() const { return -getRotation().rotate(center()); }
    void  setTranslation(Float tx, Float ty, Float tz)
    {
        auto C = -getRotation().inv().rotate(OBase {tx, ty, tz});
        memcpy(&this->x, &C.x, sizeof(Float[3]));
    }
    void setTranslation(const Point3_<Float>& tvec)
    {
        setTranslation(tvec.x, tvec.y, tvec.z);
    }
    template <typename Tp> void toMatrix(Tp rt[12]) const
    {
        auto t = getTranslation();
        for (int i = 0; i < 3; ++i) {
            (&rt[0])[i] = (Tp)((&(RBase::a1))[i]);
            (&rt[4])[i] = (Tp)((&(RBase::b1))[i]);
            (&rt[8])[i] = (Tp)((&(RBase::c1))[i]);
        }
        rt[3] = (Tp)t.x, rt[7] = (Tp)t.y, rt[11] = (Tp)t.z;
    }
    Pose_ dot(const Pose_& RT) const
    {
        return {static_cast<RBase>(*this) * static_cast<RBase>(RT),
            static_cast<RBase>(RT).inv().rotate(static_cast<OBase>(*this)) +
                static_cast<OBase>(RT)};
    }
    Pose_ inv() const
    {
        return {static_cast<RBase>(*this).inv(),
            -RBase::rotate(static_cast<OBase>(*this))};
    }
    template <typename Tp> Tp depth(const Point3_<Tp>& pt) const
    {
        return (Tp)(this->c1 * (pt.x - this->x) + this->c2 * (pt.y - this->y) +
                    this->c3 * (pt.z - this->z));
    }
    template <typename Tp>
    void depth(const Point3_<Tp>* points, Tp* depths, size_t pt_num) const
    {
        for (int i = 0; i < (int)pt_num; ++i) depths[i] = depth(points[i]);
    }
    template <typename Tp> void depth(
        const std::vector<Point3_<Tp>>& points, std::vector<Tp>& depths) const
    {
        if (points.empty()) return;
        if (points.size() != depths.size()) depths.resize(points.size());
        depth(&points[0], &depths[0], points.size());
    }
    template <typename Tp> Point3_<Tp> transform(const Point3_<Tp>& pt) const
    {
        // 考虑到精度影响的问题
        return Point3_<Tp>(RBase::rotate(
            Point3_<Float>(pt.x - OBase::x, pt.y - OBase::y, pt.z - OBase::z)));
    }
    template <typename Tp>
    void transform(const Point3_<Tp>* src, Point3_<Tp>* dst, size_t sz) const
    {
        MVS_OMP_PARALLEL_FOR
        for (int i = 0; i < (int)sz; ++i) dst[i] = transform(src[i]);
    }
#ifdef RULERMVS_USE_SSE
    void transform(const Point3f* src, Point3f* dst, size_t num) const;
#endif
    template <typename Tp>
    void transform(const std::vector<Point3_<Tp>>& points,
        std::vector<Point3_<Tp>>&                  points_t) const
    {
        if (points.empty()) return;
        if (points.size() != points_t.size()) points_t.resize(points.size());
        transform<Tp>(&points[0], &points_t[0], points.size());
    }
    /// @brief 姿态乘法运算
    Pose_& operator*=(const Pose_& rt) { return *this = dot(rt); }
    friend inline std::ostream& operator<<(std::ostream& O, const Pose_& RT)
    {
        return O << "Pose[" << static_cast<RBase>(RT) << ","
                 << static_cast<OBase>(RT) << "]";
    };
};
using Pose    = Pose_<double>;
using PoseVec = std::vector<Pose>;
template <typename T> struct CoordTraits<Pose_<T>>: CoordType<T> {};
template <typename Float> static inline Pose_<Float> operator*(
    const Pose_<Float>& rt1, const Pose_<Float>& rt2)
{
    return rt1.dot(rt2);
}
template <typename Float> static inline Pose_<Float> operator*(
    const Rotation_<Float>& R, const Pose_<Float>& rt)
{
    return {R * rt.getRotation(), rt.center()};
}
template <typename Float> static inline Pose_<Float> operator*(
    const Pose_<Float>& rt, const Rotation_<Float>& R)
{
    return {rt.getRotation() * R, R.inv().rotate(rt.center())};
}
template <typename Float>
static inline bool operator==(const Pose_<Float> RT1, const Pose_<Float>& RT2)
{
    const Float e = std::numeric_limits<Float>::epsilon();
    for (int i = 0; i < 3; ++i)
        if (std::abs((&RT1.x)[i] - (&RT2.x)[i]) > e) return false;
    return static_cast<Rotation_<Float>>(RT1) ==
           static_cast<Rotation_<Float>>(RT2);
}
template <typename Float>
static inline void toRodrigues(const Float m[9], Float rvec[3])
{
    Float a = (m[0] + m[4] + m[8] - 1) / 2;
    Float e = std::numeric_limits<Float>::epsilon();
    if (std::abs(m[1] - m[3]) < e && std::abs(m[5] - m[7]) < e &&
        std::abs(m[2] - m[6]) < e) {
        if (std::abs(m[1] + m[3]) < 0.1 && std::abs(m[5] + m[7]) < 0.1 &&
            std::abs(m[2] + m[6]) < 0.1 && a > 0.9) {
            rvec[0] = rvec[1] = rvec[2] = 0;
        } else {
            Float xx       = (m[0] + 1) / 2;
            Float yy       = (m[4] + 1) / 2;
            Float zz       = (m[8] + 1) / 2;
            Float xy       = (m[1] + m[3]) / 4;
            Float xz       = (m[2] + m[6]) / 4;
            Float yz       = (m[5] + m[7]) / 4;
            Float pi_sqrt2 = (Float)(MVS_SQRT1_2 * MVS_PI);
            if ((xx > yy) && (xx > zz)) {
                if (xx < e) {
                    rvec[0] = 0, rvec[1] = rvec[2] = pi_sqrt2;
                } else {
                    Float t = std::sqrt(xx);
                    rvec[0] = (Float)(t * MVS_PI);
                    rvec[1] = (Float)(xy / t * MVS_PI);
                    rvec[2] = (Float)(xz / t * MVS_PI);
                }
            } else if (yy > zz) {
                if (yy < e) {
                    rvec[1] = 0;
                    rvec[0] = rvec[2] = pi_sqrt2;
                } else {
                    Float t = std::sqrt(yy);
                    rvec[0] = (Float)(xy / t * MVS_PI);
                    rvec[1] = (Float)(t * MVS_PI);
                    rvec[2] = (Float)(yz / t * MVS_PI);
                }
            } else {
                if (zz < e) {
                    rvec[0] = rvec[1] = pi_sqrt2, rvec[2] = 0;
                } else {
                    Float t = std::sqrt(zz);
                    rvec[0] = (Float)(xz / t * MVS_PI);
                    rvec[1] = (Float)(yz / t * MVS_PI);
                    rvec[2] = (Float)(t * MVS_PI);
                }
            }
        }
    } else {
        a       = std::acos(a);
        Float b = (Float)(0.5 * a / std::sin(a));
        rvec[0] = b * (m[7] - m[5]);
        rvec[1] = b * (m[2] - m[6]);
        rvec[2] = b * (m[3] - m[1]);
    }
}
template <typename Float>
static inline void toRodrigues(const Rotation_<Float>& R, Float rvec[3])
{
    toRodrigues<Float>(&R.a1, rvec);
}
template <typename Float>
static inline void fromRodrigues(const Float rvec[3], Float m[9])
{
    Float a =
        std::sqrt(rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]);
    Float ct = (Float)(a == 0.0 ? 0.5 : (1.0 - std::cos(a)) / a / a);
    Float st = (Float)(a == 0.0 ? 1 : std::sin(a) / a);
    m[0]     = (Float)(1.0 - (rvec[1] * rvec[1] + rvec[2] * rvec[2]) * ct);
    m[1]     = (Float)(rvec[0] * rvec[1] * ct - rvec[2] * st);
    m[2]     = (Float)(rvec[2] * rvec[0] * ct + rvec[1] * st);
    m[3]     = (Float)(rvec[0] * rvec[1] * ct + rvec[2] * st);
    m[4]     = (Float)(1.0 - (rvec[2] * rvec[2] + rvec[0] * rvec[0]) * ct);
    m[5]     = (Float)(rvec[1] * rvec[2] * ct - rvec[0] * st);
    m[6]     = (Float)(rvec[2] * rvec[0] * ct - rvec[1] * st);
    m[7]     = (Float)(rvec[1] * rvec[2] * ct + rvec[0] * st);
    m[8]     = (Float)(1.0 - (rvec[0] * rvec[0] + rvec[1] * rvec[1]) * ct);
}
template <typename Float>
static inline decltype(auto) fromRodrigues(const Float rvec[3])
{
    Rotation_<Float> rotation;
    fromRodrigues<Float>(rvec, &rotation.a1);
    return rotation;
}
template <typename Float>
Rotation_<Float>::Rotation_(Float rx, Float ry, Float rz)
{
    Float r[3] = {rx, ry, rz};
    fromRodrigues<Float>(r, &a1);
}
template <typename Float>
static inline void toQuaternian(const Float m[9], Float quat[4])
{
    quat[0] = 1 + m[0] + m[4] + m[8];
    if (quat[0] > std::numeric_limits<Float>::epsilon()) {
        quat[0] = std::sqrt(quat[0]) / 2;
        quat[1] = (m[7] - m[5]) / (4 * quat[0]);
        quat[2] = (m[2] - m[6]) / (4 * quat[0]);
        quat[3] = (m[3] - m[1]) / (4 * quat[0]);
    } else {
        if (m[0] > m[4] && m[0] > m[8]) {
            Float s = 2 * std::sqrt(1 + m[0] - m[4] - m[8]);
            quat[1] = s / 4;
            quat[2] = (m[1] + m[3]) / s;
            quat[3] = (m[2] + m[6]) / s;
            quat[0] = (m[5] - m[7]) / s;
        } else if (m[4] > m[8]) {
            Float s = 2 * std::sqrt(1 + m[4] - m[0] - m[8]);
            quat[1] = (m[1] + m[3]) / s;
            quat[2] = s / 4;
            quat[3] = (m[5] + m[7]) / s;
            quat[0] = (m[2] - m[6]) / s;
        } else {
            Float s = 2 * std::sqrt(1 + m[8] - m[0] - m[4]);
            quat[1] = (m[2] + m[6]) / s;
            quat[2] = (m[5] + m[7]) / s;
            quat[3] = s / 4;
            quat[0] = (m[1] - m[3]) / s;
        }
    }
}
template <typename Float>
static inline void toQuaternian(const Rotation_<Float>& R, Float quat[4])
{
    toQuaternian<Float>(&R.a1, quat);
}
template <typename Float>
static inline void fromQuaternian(const Float quat[4], Float m[9])
{
    Float qq = std::sqrt(quat[0] * quat[0] + quat[1] * quat[1] +
                         quat[2] * quat[2] + quat[3] * quat[3]);
    Float qw, qx, qy, qz;
    if (qq > 0) {
        qw = quat[0] / qq;
        qx = quat[1] / qq;
        qy = quat[2] / qq;
        qz = quat[3] / qq;
    } else {
        qw = 1, qx = qy = qz = 0;
    }
    m[0] = qw * qw + qx * qx - qz * qz - qy * qy;
    m[1] = 2 * qx * qy - 2 * qz * qw;
    m[2] = 2 * qy * qw + 2 * qz * qx;
    m[3] = 2 * qx * qy + 2 * qw * qz;
    m[4] = qy * qy + qw * qw - qz * qz - qx * qx;
    m[5] = 2 * qz * qy - 2 * qx * qw;
    m[6] = 2 * qx * qz - 2 * qy * qw;
    m[7] = 2 * qy * qz + 2 * qw * qx;
    m[8] = qz * qz + qw * qw - qy * qy - qx * qx;
}
template <typename Float>
static inline decltype(auto) fromQuaternian(const Float quat[4])
{
    Rotation_<Float> rotation;
    fromQuaternian<Float>(quat, &rotation.a1);
    return rotation;
}
template <typename Float>
Rotation_<Float>::Rotation_(Float qw, Float qx, Float qy, Float qz)
{
    Float q[4] = {qw, qx, qy, qz};
    fromQuaternian<Float>(q, &a1);
}
template <typename Float>
static inline void fromEulerAxisX(Float theta, Float m[9])
{
    m[0] = 1, m[1] = 0, m[2] = 0;
    m[3] = 0, m[4] = cos(theta), m[5] = sin(theta);
    m[6] = 0, m[7] = -sin(theta), m[8] = cos(theta);
}
template <typename Float>
static inline decltype(auto) fromEulerAxisX(Float theta)
{
    Rotation_<Float> rotation;
    fromEulerAxisX<Float>(theta, &rotation.a1);
    return rotation;
}
template <typename Float>
static inline void fromEulerAxisY(Float _theta, Float _m[9])
{
    _m[0] = cos(_theta);
    _m[1] = 0;
    _m[2] = -sin(_theta);
    _m[3] = 0;
    _m[4] = 1;
    _m[5] = 0;
    _m[6] = sin(_theta);
    _m[7] = 0;
    _m[8] = cos(_theta);
}
template <typename Float>
static inline decltype(auto) fromEulerAxisY(Float _theta)
{
    Rotation_<Float> rotation;
    fromEulerAxisY<Float>(_theta, &rotation.a1);
    return rotation;
}
template <typename Float>
static inline void fromEulerAxisZ(Float _theta, Float _m[9])
{
    _m[0] = cos(_theta);
    _m[1] = sin(_theta);
    _m[2] = 0;
    _m[3] = -sin(_theta);
    _m[4] = cos(_theta);
    _m[5] = 0;
    _m[6] = 0;
    _m[7] = 0;
    _m[8] = 1;
}
template <typename Float>
static inline decltype(auto) fromEulerAxisZ(Float _theta)
{
    Rotation_<Float> rotation;
    fromEulerAxisZ<Float>(_theta, &rotation.a1);
    return rotation;
}
template <typename Float> Rotation_<Float>::Rotation_(Float theta, int axis)
{
    if (axis == 2) fromEulerAxisZ(theta, &a1);
    else if (axis == 1)
        fromEulerAxisY(theta, &a1);
    else
        fromEulerAxisX(theta, &a1);
}
template <typename Float> static inline decltype(auto) fromEulerAngles(
    const Float _angle[3], EulerAxisMode _mode = EulerAxisMode::XYZ)
{
    Rotation_<Float> rotation;
    switch (_mode) {
        case EulerAxisMode::XZY:
            rotation = fromEulerAxisY(_angle[2]);
            rotation *= fromEulerAxisZ(_angle[1]);
            rotation *= fromEulerAxisX(_angle[0]);
            break;
        case EulerAxisMode::YXZ:
            rotation = fromEulerAxisZ(_angle[2]);
            rotation *= fromEulerAxisX(_angle[1]);
            rotation *= fromEulerAxisY(_angle[0]);
            break;
        case EulerAxisMode::YZX:
            rotation = fromEulerAxisX(_angle[2]);
            rotation *= fromEulerAxisZ(_angle[1]);
            rotation *= fromEulerAxisY(_angle[0]);
            break;
        case EulerAxisMode::ZXY:
            rotation = fromEulerAxisY(_angle[2]);
            rotation *= fromEulerAxisX(_angle[1]);
            rotation *= fromEulerAxisZ(_angle[0]);
            break;
        case EulerAxisMode::ZYX:
            rotation = fromEulerAxisX(_angle[2]);
            rotation *= fromEulerAxisY(_angle[1]);
            rotation *= fromEulerAxisZ(_angle[0]);
            break;
        case EulerAxisMode::XYZ:
        default:
            rotation = fromEulerAxisZ(_angle[2]);
            rotation *= fromEulerAxisY(_angle[1]);
            rotation *= fromEulerAxisX(_angle[0]);
            break;
    }
    return rotation;
}
template <typename Float> static inline void toEulerAngles(
    const Float m[9], Float euler[3], EulerAxisMode mode = EulerAxisMode::XYZ)
{}
template <typename Float>
static inline void fromEulerAxisXYZ(const Float euler[3], Float m[9])
{
    Float alpha = euler[0], beta = euler[1], gamma = euler[2];
    m[0] = cos(gamma) * cos(beta);
    m[1] = -sin(gamma) * cos(alpha) + cos(gamma) * sin(beta) * sin(alpha);
    m[2] = sin(gamma) * sin(alpha) + cos(gamma) * sin(beta) * cos(alpha);
    m[3] = sin(gamma) * cos(beta);
    m[4] = cos(gamma) * cos(alpha) + sin(gamma) * sin(beta) * sin(alpha);
    m[5] = -cos(gamma) * sin(alpha) + sin(gamma) * sin(beta) * cos(alpha);
    m[6] = -sin(beta);
    m[7] = cos(beta) * sin(alpha);
    m[8] = cos(beta) * cos(alpha);
}
template <typename Float>
static inline decltype(auto) fromEulerAxisXYZ(const Float _angle[3])
{
    Rotation_<Float> rotation;
    fromEulerAxisXYZ<Float>(_angle, &rotation.a1);
    return rotation;
}
template <typename Float>
static inline void toEulerAxisXYZ(const Float _m[9], Float _angle[3])
{
    auto alpha   = atan2(_m[7], _m[8]);
    auto gamma   = atan2(_m[3], _m[0]);
    auto cosbeta = (_m[0] * cos(gamma) + _m[3] * sin(gamma) +
                       _m[7] * sin(alpha) + _m[8] * cos(alpha)) /
                   2;
    auto beta = atan2(-_m[6], cosbeta);
    _angle[0] = alpha;
    _angle[1] = beta;
    _angle[2] = gamma;
}
template <typename Float> static inline void toEulerAxisXYZ(
    const Rotation_<Float>& rotation, Float angle[3])
{
    toEulerAxisXYZ(&rotation.a1, angle);
}
template <typename Float>
static inline void updatePoseFromX(Pose_<Float>& rt, const Float X[6])
{
    Rotation_<Float> R = fromEulerAxisXYZ(&X[0]);
    rt *= Pose_<Float>(R, -R.inv().rotate(Point3_<Float> {X[3], X[4], X[5]}));
}

template <typename Tp, typename Float>
static inline void rotatePoints(const Rotation_<Float>& r,
    const Point3_<Tp>* src, Point3_<Tp>* dst, size_t num)
{
    if (src && dst && num) r.rotate(src, dst, num);
}

template <typename Tp, typename Float>
static inline void rotatePoints(const Rotation_<Float>& r,
    const std::vector<Point3_<Tp>>& src, std::vector<Point3_<Tp>>& dst)
{
    if (src.empty()) return;
    if (dst.size() != src.size()) dst.resize(src.size());
    rotatePoints(r, src.data(), dst.data(), src.size());
}

template <typename Tp, typename Float> static inline void transformPoints(
    const Pose_<Float>& p, const Point3_<Tp>* src, Point3_<Tp>* dst, size_t num)
{
    if (src && dst && num) p.transform(src, dst, num);
}

template <typename Tp, typename Float>
static inline void transformPoints(const Pose_<Float>& p,
    const std::vector<Point3_<Tp>>& src, std::vector<Point3_<Tp>>& dst)
{
    if (src.empty()) return;
    if (dst.size() != src.size()) dst.resize(src.size());
    transformPoints(p, src.data(), dst.data(), src.size());
}

// RULERMVS_JSON_IO_EXPORT(Pose);
// RULERMVS_JSON_IO_EXPORT(PoseVec);
// RULERMVS_JSON_IO_EXPORT(Rotation);

// /// @brief R*[x,y,z]'的指令集加速
// template <typename Float> void Rotation_<Float>::rotate(
//     const Point3_<float>* src, Point3_<float>* dst, size_t num) const
// {
//     // 如果平台不支持指令集加速，则采用一般方法.
//     if (!checkHardwareSupport(SIMDMode::MVS_SSE3)) {
//         this->rotate<float>(src, dst, num);
//         return;
//     }
//     __m128 R[9];
//     for (int i = 0; i < 9; ++i) R[i] = _mm_set1_ps((float)(&(this->a1))[i]);

//     MVS_OMP_PARALLEL_FOR
//     for (int i = 0; i < (int)num >> 2; ++i) {
//         auto*  src_ptr = &src[i << 2];
//         auto*  dst_ptr = &dst[i << 2];
//         __m128 X =
//             _mm_set_ps(src_ptr[3].x, src_ptr[2].x, src_ptr[1].x,
//             src_ptr[0].x);
//         __m128 Y =
//             _mm_set_ps(src_ptr[3].y, src_ptr[2].y, src_ptr[1].y,
//             src_ptr[0].y);
//         __m128 Z =
//             _mm_set_ps(src_ptr[3].z, src_ptr[2].z, src_ptr[1].z,
//             src_ptr[0].z);
//         __m128 Xt =
//             _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, R[0]), _mm_mul_ps(Y, R[1])),
//                 _mm_mul_ps(Z, R[2]));
//         __m128 Yt =
//             _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, R[3]), _mm_mul_ps(Y, R[4])),
//                 _mm_mul_ps(Z, R[5]));
//         __m128 Zt =
//             _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, R[6]), _mm_mul_ps(Y, R[7])),
//                 _mm_mul_ps(Z, R[8]));
//         for (int j = 0; j < 4; ++j) {
//             dst_ptr[j].x = ((float*)&Xt)[j];
//             dst_ptr[j].y = ((float*)&Yt)[j];
//             dst_ptr[j].z = ((float*)&Zt)[j];
//         }
//     }
//     for (size_t i = (num >> 2) << 2; i < num; ++i) dst[i] = rotate(src[i]);
// }

// /// @brief RT*[x,y,z]'的SSE加速
// template <typename Float>
// void Pose_<Float>::transform(const Point3f* src, Point3f* dst, size_t num)
// const
// {
//     // 如果平台不支持指令集加速，则采用一般方法.
//     if (!checkHardwareSupport(SIMDMode::MVS_SSE3)) {
//         this->transform<float>(src, dst, num);
//         return;
//     }
//     __m128 RT[12];
//     for (int i = 0; i < 3; ++i) {
//         (&RT[0])[i] = _mm_set1_ps((float)((&(this->a1))[i]));
//         (&RT[4])[i] = _mm_set1_ps((float)((&(this->b1))[i]));
//         (&RT[8])[i] = _mm_set1_ps((float)((&(this->c1))[i]));
//     }
//     auto t = this->getTranslation();
//     RT[3]  = _mm_set1_ps((float)t.x);
//     RT[7]  = _mm_set1_ps((float)t.y);
//     RT[11] = _mm_set1_ps((float)t.z);

//     MVS_OMP_PARALLEL_FOR
//     for (int i = 0; i < (int)num >> 2; ++i) {
//         auto*  src_ptr = &src[i << 2];
//         auto*  dst_ptr = &dst[i << 2];
//         __m128 X =
//             _mm_set_ps(src_ptr[3].x, src_ptr[2].x, src_ptr[1].x,
//             src_ptr[0].x);
//         __m128 Y =
//             _mm_set_ps(src_ptr[3].y, src_ptr[2].y, src_ptr[1].y,
//             src_ptr[0].y);
//         __m128 Z =
//             _mm_set_ps(src_ptr[3].z, src_ptr[2].z, src_ptr[1].z,
//             src_ptr[0].z);
//         __m128 Xt = _mm_add_ps(
//             _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, RT[0]), _mm_mul_ps(Y,
//             RT[1])),
// _mm_mul_ps(Z, RT[2])),
//             RT[3]);
//         __m128 Yt = _mm_add_ps(
//             _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, RT[4]), _mm_mul_ps(Y,
//             RT[5])),
//                 _mm_mul_ps(Z, RT[6])),
//             RT[7]);
//         __m128 Zt = _mm_add_ps(
//             _mm_add_ps(_mm_add_ps(_mm_mul_ps(X, RT[8]), _mm_mul_ps(Y,
//             RT[9])),
//                 _mm_mul_ps(Z, RT[10])),
//             RT[11]);
//         for (int j = 0; j < 4; ++j) {
//             dst_ptr[j].x = ((float*)&Xt)[j];
//             dst_ptr[j].y = ((float*)&Yt)[j];
//             dst_ptr[j].z = ((float*)&Zt)[j];
//         }
//     }
//     for (size_t i = (num >> 2) << 2; i < num; ++i) dst[i] =
//     transform(src[i]);
// }

/// @brief 畸变模型的枚举类型
enum class DistorType : int { None, Radial, Brown };
/// @brief 占位，对应无相机畸变
/// @tparam Float 浮点类型
/// @tparam D 畸变类型
template <DistorType D = DistorType::None, typename Float = double>
struct DistorModel_ {
    DistorModel_() {}
    DistorModel_(Float...) {}
    DistorModel_(const DistorModel_&) {}
    template <typename Tp,
        typename = std::enable_if_t<!std::is_same<Float, Tp>::value>>
    explicit DistorModel_(const DistorModel_<D, Tp>&)
    {}
    DistorModel_& operator=(const DistorModel_&) { return *this; }
    template <typename Tp>
    const Point3_<Tp>& distor(const Point3_<Tp>& pt) const noexcept
    {
        return static_cast<const Point3_<Tp>&>(pt);
    }
    template <typename Tp>
    const Point3_<Tp>& undistor(const Point3_<Tp>& pt) const noexcept
    {
        return static_cast<const Point3_<Tp>&>(pt);
    }
    friend inline std::ostream& operator<<(std::ostream& o, const DistorModel_&)
    {
        return o << "DistorNone";
    };
    const Float k1 = 0;  ///< 占位符
};
template <typename Float> struct DistorModel_<DistorType::Radial, Float> {
    DistorModel_() : k1(0), k2(0), k3(0) {}
    DistorModel_(const DistorModel_<DistorType::None, Float>& d)
        : k1(d.k1), k2(0), k3(0)
    {}
    DistorModel_(const DistorModel_& d) : k1(d.k1), k2(d.k2), k3(d.k3) {}
    template <typename Tp,
        typename = std::enable_if_t<!std::is_same<Float, Tp>::value>>
    explicit DistorModel_(const DistorModel_<DistorType::Radial, Tp>& d)
        : k1((Float)d.k1), k2((Float)d.k2), k3((Float)d.k3)
    {}
    DistorModel_(Float k1, Float k2 = 0, Float k3 = 0) : k1(k1), k2(k2), k3(k3)
    {}
    template <typename Tp> Point3_<Tp> distor(const Point3_<Tp>& pt) const
    {
        if (!pt.z) return {0, 0, 0};
        Float x  = pt.x / pt.z;
        Float y  = pt.y / pt.z;
        Float r  = x * x + y * y;
        Float rr = r * r;
        Float kd = 1 + k1 * r + k2 * rr + k3 * rr * r;
        return {(Tp)(x * kd), (Tp)(y * kd), (Tp)1};
    }
    template <typename Tp>
    Point3_<Tp> undistor(const Point3_<Tp>& pt, int max_iter = 5) const
    {
        if (!pt.z) return {0, 0, 0};
        Float x0 = pt.x / pt.z;
        Float y0 = pt.y / pt.z;
        Float x = x0, y = y0;
        for (int i = 0; i < max_iter; ++i) {
            Float r  = x * x + y * y;
            Float rr = r * r;
            Float kd = 1 + k1 * r + k2 * rr + k3 * rr * r;

            Float b1 = kd * x - x0;
            Float b2 = kd * y - y0;
            Float J11 =
                kd + x * (2 * k1 * x + 4 * k2 * x * r + 6 * k3 * x * rr);
            Float J12 = x * (2 * k1 * y + 4 * k2 * y * r + 6 * k3 * y * rr);
            Float J21 = y * (2 * k1 * x + 4 * k2 * x * r + 6 * k3 * x * rr);
            Float J22 =
                kd + y * (2 * k1 * y + 4 * k2 * y * r + 6 * k3 * y * rr);

            Float A11 = J11 * J11 + J21 * J21;
            Float A12 = J11 * J12 + J21 * J22;
            Float A21 = J12 * J11 + J22 * J21;
            Float A22 = J12 * J12 + J22 * J22;
            Float B1  = J11 * b1 + J12 * b2;
            Float B2  = J21 * b1 + J22 * b2;
            Float det = A11 * A22 - A12 * A21;
            Float X1  = (A22 * B1 - A12 * B2) / det;
            Float X2  = (A11 * B2 - A21 * B1) / det;

            x -= X1;
            y -= X2;
        }
        return {(Tp)x, (Tp)y, (Tp)1};
    }
    friend inline std::ostream& operator<<(
        std::ostream& o, const DistorModel_& d)
    {
        return o << "DistorRadial(" << d.k1 << "," << d.k2 << "," << d.k3
                 << ")";
    };

    Float k1;  ///< 径向畸变K1
    Float k2;  ///< 径向畸变K2
    Float k3;  ///< 径向畸变K3
};
template <typename Float> struct DistorModel_<DistorType::Brown, Float> {
    DistorModel_() : k1(0), k2(0), p1(0), p2(0), k3(0) {}
    DistorModel_(const DistorModel_& d)
        : k1(d.k1), k2(d.k2), p1(d.p1), p2(d.p2), k3(d.k3)
    {}
    template <typename Tp,
        typename = std::enable_if_t<!std::is_same<Float, Tp>::value>>
    explicit DistorModel_(const DistorModel_<DistorType::Brown, Tp>& d)
        : k1((Float)d.k1)
        , k2((Float)d.k2)
        , p1((Float)d.p1)
        , p2((Float)d.p2)
        , k3((Float)d.k3)
    {}
    DistorModel_(
        Float k1, Float k2 = 0, Float p1 = 0, Float p2 = 0, Float k3 = 0)
        : k1(k1), k2(k2), p1(p1), p2(p2), k3(k3)
    {}
    template <typename Tp> Point3_<Tp> distor(const Point3_<Tp>& pt) const
    {
        Float x   = pt.x / pt.z;
        Float y   = pt.y / pt.z;
        Float xx  = x * x;
        Float yy  = y * y;
        Float xy  = x * y;
        Float r   = xx + yy;
        Float rr  = r * r;
        Float kd  = 1 + k1 * r + k2 * rr + k3 * rr * r;
        Float dox = 2 * p1 * xy + p2 * (r + 2 * xx);
        Float doy = p1 * (r + 2 * yy) + 2 * p2 * xy;
        return {(Tp)(x * kd + dox), (Tp)(y * kd + doy), (Tp)1};
    }
    template <typename Tp>
    Point3_<Tp> undistor(const Point3_<Tp>& pt, int max_iter = 5) const
    {
        Float x0 = pt.x / pt.z;
        Float y0 = pt.y / pt.z;
        Float x = x0, y = y0;
        for (int i = 0; i < max_iter; ++i) {
            Float xx  = x * x;
            Float xy  = x * y;
            Float yy  = y * y;
            Float r   = xx + yy;
            Float rr  = r * r;
            Float kd  = 1 + k1 * r + k2 * rr + k3 * rr * r;
            Float dox = 2 * p1 * xy + p2 * (r + 2 * xx);
            Float doy = p1 * (r + 2 * yy) + 2 * p2 * xy;

            Float b1  = x * kd + dox - x0;
            Float b2  = y * kd + doy - y0;
            Float J11 = kd +
                        x * (2 * k1 * x + 4 * k2 * x * r + 6 * k3 * x * rr) +
                        2 * p1 * y + 6 * p2 * x;
            Float J12 = x * (2 * k1 * y + 4 * k2 * y * r + 6 * k3 * y * rr) +
                        2 * p1 * x + 2 * p2 * y;
            Float J21 = y * (2 * k1 * x + 4 * k2 * x * r + 6 * k3 * x * rr) +
                        2 * p1 * x + 2 * p2 * y;
            Float J22 = kd +
                        y * (2 * k1 * y + 4 * k2 * y * r + 6 * k3 * y * rr) +
                        6 * p1 * y + 2 * p2 * x;

            Float A11 = J11 * J11 + J21 * J21;
            Float A12 = J11 * J12 + J21 * J22;
            Float A21 = J12 * J11 + J22 * J21;
            Float A22 = J12 * J12 + J22 * J22;
            Float B1  = J11 * b1 + J12 * b2;
            Float B2  = J21 * b1 + J22 * b2;
            Float det = A11 * A22 - A12 * A21;
            Float X1  = (A22 * B1 - A12 * B2) / det;
            Float X2  = (A11 * B2 - A21 * B1) / det;

            x -= X1;
            y -= X2;
        }
        return {(Tp)x, (Tp)y, (Tp)1};
    }
    friend inline std::ostream& operator<<(
        std::ostream& o, const DistorModel_& d)
    {
        return o << "DistorBrown(" << d.k1 << "," << d.k2 << "," << d.p1 << ","
                 << d.p2 << "," << d.k3 << ")";
    };

    Float k1;  ///< 经向畸变K1
    Float k2;  ///< 经向畸变K2
    Float p1;  ///< 切向畸变P1
    Float p2;  ///< 切向畸变P2
    Float k3;  ///< 经向畸变K3
};
typedef DistorModel_<DistorType::None, double>   NoneDistor;
typedef DistorModel_<DistorType::Brown, double>  BrownDistor;
typedef DistorModel_<DistorType::Radial, double> RadialDistor;

/// @brief 相机类型枚举
enum class CameraType : int {
    SixBox,      ///< 六面体
    Fisheye,     ///< 鱼眼
    Pinhole,     ///< 针孔相机模型
    Spherical,   ///< 全景
    SkewPinhole  ///< 带斜率针孔相机模型
};

/// @brief 相机模型
template <CameraType C = CameraType::Pinhole, typename Float = double>
struct CameraModel_ {
    CameraModel_() : fx(0), fy(0), cx(0), cy(0) {}
    CameraModel_(Float f, Float cx, Float cy) : fx(f), fy(f), cx(cx), cy(cy) {}
    CameraModel_(Float fx, Float fy, Float cx, Float cy)
        : fx(fx), fy(fy), cx(cx), cy(cy)
    {}
    CameraModel_(const CameraModel_& p) : fx(p.fx), fy(p.fy), cx(p.cx), cy(p.cy)
    {}
    template <typename Tp,
        typename = std::enable_if_t<!std::is_same<Float, Tp>::value>>
    explicit CameraModel_(const CameraModel_<CameraType::Pinhole, Tp>& p)
        : fx((Float)p.fx), fy((Float)p.fy), cx((Float)p.cx), cy((Float)p.cy)
    {}
    template <typename Tp> Point2_<Tp> point2uv(const Point3_<Tp>& pt) const
    {
        return {(Tp)(fx * pt.x / pt.z + cx), (Tp)(fy * pt.y / pt.z + cy)};
    }
    template <typename Tp> Point3_<Tp> uv2point(const Point2_<Tp>& uv) const
    {
        return {(Tp)((uv.x - cx) / fx), (Tp)((uv.y - cy) / fy), (Tp)1};
    }
    // friend CameraModel_ operator/(
    //     const CameraModel_& cam, coord_traits_t<Float> s)
    // {
    //     return {cam.fx / s, cam.fy / s, cam.cx / s, cam.cy / s};
    // }
    // friend CameraModel_ operator*(
    //     const CameraModel_& cam, const Scalar2_<Float>& s)
    // {
    //     return {cam.fx * s[0], cam.fy * s[1], cam.cx * s[0], cam.cy * s[1]};
    // }
    // friend CameraModel_ operator/(
    //     const CameraModel_& cam, const Scalar2_<Float>& s)
    // {
    //     return {cam.fx / s[0], cam.fy / s[1], cam.cx / s[0], cam.cy / s[1]};
    // }
    friend inline std::ostream& operator<<(
        std::ostream& O, const CameraModel_& cam)
    {
        return O << "Pinhole(" << cam.fx << "," << cam.fy << "," << cam.cx
                 << "," << cam.cy << ")";
    };

    Float fx;  ///< X方向的焦距
    Float fy;  ///< Y方向的焦距
    Float cx;  ///< X方向的像主点坐标
    Float cy;  ///< Y方向的像主点坐标
};

/// @brief 带斜率针孔相机模型
template <typename Float> struct CameraModel_<CameraType::SkewPinhole, Float> {
    CameraModel_() : fx(0), fy(0), cx(0), cy(0), sk(0) {}
    CameraModel_(Float fx, Float fy, Float cx, Float cy, Float sk = 0)
        : fx(fx), fy(fy), cx(cx), cy(cy), sk(sk)
    {}
    CameraModel_(const CameraModel_& p)
        : fx(p.fx), fy(p.fy), cx(p.cx), cy(p.cy), sk(p.sk)
    {}
    template <typename Tp,
        typename = std::enable_if_t<!std::is_same<Float, Tp>::value>>
    explicit CameraModel_(const CameraModel_<CameraType::Pinhole, Tp>& p)
        : fx((Float)p.fx)
        , fy((Float)p.fy)
        , cx((Float)p.cx)
        , cy((Float)p.cy)
        , sk((Float)p.sk)
    {}
    template <typename Tp> Point2_<Tp> point2uv(const Point3_<Tp>& pt) const
    {
        return {(Tp)(fx * pt.x / pt.z + sk * pt.y / pt.z + cx),
            (Tp)(fy * pt.y / pt.z + cy)};
    }
    template <typename Tp> Point3_<Tp> uv2point(const Point2_<Tp>& uv) const
    {
        double v = (uv.y - cy) / fy;
        double u = (uv.x - sk * v - cx) / fx;
        return {(Tp)u, (Tp)v, (Tp)1};
    }

    // friend CameraModel_ operator*(
    //     const CameraModel_& cam, coord_traits_t<Float> s)
    // {
    //     return {cam.fx * s, cam.fy * s, cam.cx * s, cam.cy * s};
    // }
    // friend CameraModel_ operator/(
    //     const CameraModel_& cam, coord_traits_t<Float> s)
    // {
    //     return {cam.fx / s, cam.fy / s, cam.cx / s, cam.cy / s};
    // }
    // friend CameraModel_ operator*(
    //     const CameraModel_& cam, const Scalar2_<Float>& s)
    // {
    //     return {cam.fx * s[0], cam.fy * s[1], cam.cx * s[0], cam.cy * s[1],
    //         cam.sk * s[0]};
    // }
    // friend CameraModel_ operator/(
    //     const CameraModel_& cam, const Scalar2_<Float>& s)
    // {
    //     assert(s[0] != 0 && s[1] != 0);
    //     return {cam.fx / s[0], cam.fy / s[1], cam.cx / s[0], cam.cy / s[1],
    //         cam.sk / s[0]};
    // }
    friend inline std::ostream& operator<<(
        std::ostream& O, const CameraModel_& cam)
    {
        return O << "SkewPinhole(" << cam.fx << "," << cam.fy << "," << cam.cx
                 << "," << cam.cy << "," << cam.sk << ")";
    };

    Float fx;  ///< X方向的焦距
    Float fy;  ///< Y方向的焦距
    Float cx;  ///< X方向的像主点坐标
    Float cy;  ///< Y方向的像主点坐标
    Float sk;  ///< 斜率，微距相机中存在,或者标定板不正或提点精度不高容易造成.
};

template <typename Float> struct CameraModel_<CameraType::SixBox, Float> {
    Point2_<Float> point2uv(const Point3_<Float>& pt) const { return {0, 0}; }
    Point3_<Float> uv2point(const Point2_<Float>& uv) const
    {
        return {0, 0, 0};
    }

    Float square;
};
template <typename Float> struct CameraModel_<CameraType::Spherical, Float> {
    CameraModel_() {}
    Point2_<Float> point2uv(const Point3_<Float>& pt) const { return {0, 0}; }
    Point3_<Float> uv2point(const Point2_<Float>& uv) const
    {
        return {0, 0, 0};
    }

    Float ppx;
    Float ppy;
};
typedef CameraModel_<CameraType::SixBox, double>      SixBoxModel;
typedef CameraModel_<CameraType::Pinhole, double>     PinholeModel;
typedef CameraModel_<CameraType::Spherical, double>   SphericalModel;
typedef CameraModel_<CameraType::SkewPinhole, double> SkewPinholeModel;

template <typename Float, CameraType C, DistorType D = DistorType::None>
struct Camera_: CameraModel_<C, Float>, DistorModel_<D, Float> {
    static_assert(std::is_floating_point<Float>::value,
        "Struct Camera_ Can Only Be Instantiated With Float Point Type.");
    typedef CameraModel_<C, Float> CameraBase;
    typedef DistorModel_<D, Float> DistorBase;
    Camera_() : CameraBase(), DistorBase(), width(0), height(0) {}
    // explicit Camera_(Size sz)
    //     : CameraBase(), DistorBase(), width(sz.width), height(sz.height)
    // {}
    Camera_(const CameraBase& c) : CameraBase(c), DistorBase() {}
    Camera_(const CameraBase& c, const DistorBase& d)
        : CameraBase(c), DistorBase(d)
    {}
    Camera_(const Camera_& cam)
        : CameraBase(cam), DistorBase(cam), width(cam.width), height(cam.height)
    {}
    Camera_(const std::string& path) : Camera_() { load(path); }
    template <typename Tp,
        typename = std::enable_if_t<!IsSame<Float, Tp>::value>>
    explicit Camera_(const Camera_<Tp, C, D>& cam)
        : width(cam.width), height(cam.height), CameraBase(cam), DistorBase(cam)
    {}
    bool load(const std::string& path);
    bool save(const std::string& path) const;
    // Size size() const { return {width, height}; }
    /// @brief 重置相机对应影像分辨率的参数.
    /// @return 返回自身引用
    // Camera_& resize(Size sz);
    /// @brief 返回无畸变内参
    Camera_<Float, C, DistorType::None> nodistor() const
    {
        return {this->size(), static_cast<CameraBase>(*this)};
    }
    template <typename = typename EnableIf<C == CameraType::SkewPinhole>::type>
    Camera_<Float, CameraType::Pinhole, D> noskew() const
    {
        return {this->size(),
            {CameraBase::fx, CameraBase::fy, CameraBase::cx, CameraBase::cy},
            static_cast<DistorBase>(*this)};
    }
    template <typename Tp> void uv2points(
        const Point2_<Tp>* uvs, Point3_<Tp>* pts, size_t pt_num) const
    {
        assert(uvs != nullptr && pts != nullptr && pt_num > 0);
        for (size_t i = 0; i < pt_num; ++i)
            pts[i] = CameraBase::uv2point(uvs[i]);
    }
    template <typename Tp> void uv2points(const std::vector<Point2_<Tp>>& uvs,
        std::vector<Point3_<Tp>>& pts) const
    {
        if (uvs.empty()) return;
        if (uvs.size() != pts.size()) pts.resize(uvs.size());
        uv2points(&uvs[0], &pts[0], uvs.size());
    }
    template <typename Tp> void point2uvs(
        const Point3_<Tp>* pts, Point2_<Tp>* uvs, size_t pt_num) const
    {
        assert(pts != nullptr && uvs != nullptr && pt_num > 0);
        for (size_t i = 0; i < pt_num; ++i)
            uvs[i] = CameraBase::point2uv(pts[i]);
    }
    template <typename Tp> void point2uvs(const std::vector<Point3_<Tp>>& pts,
        std::vector<Point2_<Tp>>& uvs) const
    {
        if (pts.empty()) return;
        if (pts.size() != uvs.size()) uvs.resize(pts.size());
        point2uvs(&pts[0], &uvs[0], pts.size());
    }
    template <typename Tp> Point3_<Tp> project(const Point2_<Tp>& uv) const
    {
        return DistorBase::undistor(CameraBase::uv2point(uv));
    }
    template <typename Tp>
    void project(const Point2_<Tp>* uvs, Point3_<Tp>* pts, size_t pt_num) const
    {
        assert(uvs != nullptr && pts != nullptr && pt_num > 0);
        for (size_t i = 0; i < pt_num; ++i) pts[i] = project(uvs[i]);
    }
    template <typename Tp> void project(const std::vector<Point2_<Tp>>& uvs,
        std::vector<Point3_<Tp>>& pts) const
    {
        if (uvs.empty()) return;
        if (uvs.size() != pts.size()) pts.resize(uvs.size());
        project(&uvs[0], &pts[0], uvs.size());
    }
    template <typename Tp> Point2_<Tp> reproject(const Point3_<Tp>& pt) const
    {
        return CameraBase::point2uv(DistorBase::distor(pt));
    }
    template <typename Tp> void reproject(
        const Point3_<Tp>* pts, Point2_<Tp>* uvs, size_t pt_num) const
    {
        assert(pts != nullptr && uvs != nullptr && pt_num > 0);
        for (size_t i = 0; i < pt_num; ++i) uvs[i] = reproject(pts[i]);
    }
    template <typename Tp> void reproject(const std::vector<Point3_<Tp>>& pts,
        std::vector<Point2_<Tp>>& uvs) const
    {
        if (pts.empty()) return;
        if (uvs.size() != pts.size()) uvs.resize(pts.size());
        reproject(&pts[0], &uvs[0], uvs.size());
    }
    friend inline std::ostream& operator<<(std::ostream& o, const Camera_& cam)
    {
        return o << "Camera[ " << static_cast<CameraBase>(cam) << " "
                 << static_cast<DistorBase>(cam) << " " << cam.size() << " ]";
    };

    // int width, height;
};
template <typename Float, DistorType D> using CameraP_ =
    Camera_<Float, CameraType::Pinhole, D>;
template <typename Float, DistorType D> using CameraSkewP_ =
    Camera_<Float, CameraType::SkewPinhole, D>;
typedef CameraP_<double, DistorType::None>       CameraP;
typedef CameraP_<double, DistorType::Brown>      CameraPB;
typedef CameraP_<double, DistorType::Radial>     CameraPR;
typedef CameraSkewP_<double, DistorType::None>   CameraSkewP;
typedef CameraSkewP_<double, DistorType::Brown>  CameraSkewPB;
typedef CameraSkewP_<double, DistorType::Radial> CameraSkewPR;

// template <typename Float, CameraType C, DistorType D>
// static inline decltype(auto) resizeCamera(
//     const Camera_<Float, C, D>& cam, Size sz)
// {
//     using CameraT = Camera_<Float, C, D>;
//     if (cam.size() == sz) return CameraT {cam};
//     Scalar2_<Float> s = {1.0f, 1.0f};
//     if (cam.width > 0 && cam.height > 0) {
//         s[0] = static_cast<Float>((double)sz.width / cam.width);
//         s[1] = static_cast<Float>((double)sz.height / cam.height);
//     }
//     return CameraT {sz, static_cast<typename CameraT::CameraBase>(cam) * s,
//         static_cast<typename CameraT::DistorBase>(cam)};
// }

// template <typename Float, CameraType C, DistorType D>
// Camera_<Float, C, D>& Camera_<Float, C, D>::resize(Size sz)
// {
//     return *this = resizeCamera(*this, sz);
// }
template <typename Float, CameraType C, DistorType D>
static inline Camera_<Float, C, D> operator*(
    const Camera_<Float, C, D>& cam, double s)
{
    return resizeCamera(cam, Size_<int> {cam.width * s, cam.height * s});
}
// static inline Camera_<Float, C, D> operator*(
//     const Camera_<Float, C, D>& cam, coord_traits_t<Float> s)
// {
//     return resizeCamera(cam, Size_<int> {static_cast<int>(cam.width * s),
//                                  static_cast<int>(cam.height * s)});
// }
template <typename Float, CameraType C, DistorType D>
static inline Camera_<Float, C, D> operator/(
    const Camera_<Float, C, D>& cam, double s)
{
    return resizeCamera(cam, Size_<int> {cam.width / s, cam.height / s});
}

// 畸变矫正
template <typename Tp, typename Float, CameraType C, DistorType D>
void undistorPoints(const Camera_<Float, C, D>& cam, const Point2_<Tp>* uvs,
    Point2_<Tp>* un_uvs, size_t sz)
{
    for (size_t i = 0; i < sz; ++i)
        un_uvs[i] = cam.point2uv(cam.undistor(cam.uv2point(uvs[i])));
}
template <typename Tp, typename Float, CameraType C, DistorType D>
void undistorPoints(const Camera_<Float, C, D>& cam,
    const std::vector<Point2_<Tp>>& uvs, std::vector<Point2_<Tp>>& un_uvs)
{
    if (uvs.empty()) return;
    if (uvs.size() != un_uvs.size()) un_uvs.resize(uvs.size());
    undistorPoints<Tp, Float, C, D>(cam, &uvs[0], &un_uvs[0], uvs.size());
}
template <typename Tp, typename Float, CameraType C, DistorType D>
decltype(auto) undistorPoints(
    const Camera_<Float, C, D>& cam, const std::vector<Point2_<Tp>>& uvs)
{
    typename std::remove_cv_t<std::remove_reference_t<decltype(uvs)>> un_uvs;
    undistorPoints<Tp, Float, C, D>(cam, uvs, un_uvs);
    return un_uvs;
}

// 从JSON文件读取或保存相机内参
// RULERMVS_JSON_IO_EXPORT(CameraP);
// RULERMVS_JSON_IO_EXPORT(CameraPB);

// template <typename Float, CameraType C, DistorType D>
// bool Camera_<Float, C, D>::load(ConstStr& path)
// {
//     return JSONReader_<Camera_<Float, C, D>>::read(path, *this);
// }
// template <typename Float, CameraType C, DistorType D>
// bool Camera_<Float, C, D>::save(ConstStr& path) const
// {
//     return JSONWriter_<Camera_<Float, C, D>>::write(path, *this);
// }
// template <typename Float, CameraType C, DistorType D>
// struct JSONReader_<Camera_<Float, C, D>,
//     std::enable_if_t<has_json_reader<Camera_<Float, C, D>>::value>> {
//     static inline bool read(ConstStr& path, Camera_<Float, C, D>& camera)
//     {
//         return readJson(path, camera);
//     }
// };
// template <typename Float, CameraType C, DistorType D>
// struct JSONWriter_<Camera_<Float, C, D>,
//     std::enable_if_t<has_json_writer<Camera_<Float, C, D>>::value>> {
//     static inline bool write(ConstStr& path, const Camera_<Float, C, D>&
//     camera)
//     {
//         return writeJson(path, camera);
//     }
// };

}  // namespace minicam

#endif  // _MINICAM_MINICAM_HPP_
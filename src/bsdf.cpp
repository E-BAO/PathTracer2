#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {
  Vector3D z = Vector3D(n.x, n.y, n.z);
  Vector3D h = z;
  if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
  else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
  else h.z = 1.0;

  z.normalize();
  Vector3D y = cross(h, z);
  y.normalize();
  Vector3D x = cross(z, y);
  x.normalize();

  o2w[0] = x;
  o2w[1] = y;
  o2w[2] = z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 1.2
  // Using BSDF::reflect(), implement sample_f for a mirror surface
    *pdf = 1.0f;
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
  // TODO: 2.2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
    float tan_theta_h = sin_theta(h) / cos_theta(h);
    return exp(- pow(tan_theta_h, 2)/pow(alpha, 2)) / (PI * pow(alpha, 2) * pow(cos_theta(h), 4));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
  // TODO: 2.3
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Spectrum.
    Spectrum eta2_k2 = eta * eta + k * k;
    Spectrum Rs = ((eta2_k2) - 2. * eta * cos_theta(wi) + pow(cos_theta(wi), 2)) \
    / ((eta2_k2) + 2. * eta * cos_theta(wi) + pow(cos_theta(wi), 2));
    Spectrum Rp = ((eta2_k2) * pow(cos_theta(wi), 2) - 2. * eta * cos_theta(wi) + 1.) \
    / ((eta2_k2) * pow(cos_theta(wi), 2) + 2. * eta * cos_theta(wi) + 1.);
    
    return (Rs + Rp) * 0.5;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO: 2.1
  // Implement microfacet model here
    if(wi.z <= 0. || wo.z <= 0.){
        return Spectrum();
    }
    Vector3D h = (wo + wi);
    h.normalize();
    return (F(wi) * G(wo, wi) * D(h)) / 4. / wo.z / wi.z;
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 2.4
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.
    
//    *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
//    return MicrofacetBSDF::f(wo, *wi);
    
    Vector2D r = sampler.get_sample();
    
    float tan_theta_h = sqrt(-pow(alpha, 2) * log(1. - r[0]));
    float theta_h = atan(tan_theta_h);
    float phi_h = 2. * PI * r[1];
    Vector3D h(cos(phi_h) * sin(theta_h), sin(phi_h) * sin(theta_h), cos(theta_h));
    *wi = -wo + 2. * dot(wo, h) * h;
    
    if(wi->z <= 0. || wo.z <= 0.){
        *pdf = 0.;
        return Spectrum();
    }
    
    float p_h = 2. * sin(theta_h)* exp(- pow(tan_theta_h, 2)/ pow(alpha, 2)) / pow(cos(theta_h), 3) / pow(alpha, 2);
    
    float p_phi = 0.5 / PI;
    float p_w_h = p_phi * p_h / sin(theta_h);
    float p_w_wi = p_w_h / 4.0 / dot(*wi, h);
    
    *pdf = p_w_wi;//min(p_w_wi, 1.0f - EPS_F);
    return f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO: 1.4
  // Compute Fresnel coefficient and either reflect or refract based on it.

    if(!refract(wo, wi, ior)){
        reflect(wo, wi);
        *pdf = 1.0;
        return reflectance / abs_cos_theta(*wi);
    }else{
        float R0 = pow((1. - ior) / (1. + ior),2);
        float R = R0 + (1. - R0)*pow(1 - abs_cos_theta(wo), 5);
        if(coin_flip(R)){
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }else{
            *pdf = 1.0 - R;
            float eta = wo.z < 0? ior: 1. / ior;
            return (1.-R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
        }
    }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO: 1.1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO: 1.3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    float eta = wo.z < 0? ior: 1. / ior;
    float costheta = wo.z;
    float sq = 1.0 - eta * eta * (1. - costheta * costheta);
    if(sq < 0)
        return false;
    else{
        *wi = Vector3D(-eta * wo.x, -eta * wo.y, wo.z >= 0? - sqrt(sq): sqrt(sq));
        return true;
    }
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL

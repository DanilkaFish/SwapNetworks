import { errorResponce } from "../util/utils";
import * as userServices from "../service/user.service";
import { NextFunction } from "express";

export const register = async (req: any, res: any) => {
    const { email, password } = req.body;
    try {
        const userInfo = await userServices.createUser(email, password);
        res.status(200).send(userInfo);
    } catch (error) {
        errorResponce(error, res);
    }
};

export const login = async (req: any, res: any) => {
    try {
        const { email, password } = req.body;
        const userInfo = await userServices.grantAccess(email, password);
        res.status(200).send(userInfo);
    } catch (error) {
        errorResponce(error, res);
    }
};

export const refreshToken = async (req: any, res: any) => {
    try {
        const { email, rt } = req.body;
        const userInfo = await userServices.refreshToken(email, rt);
        res.status(200).send(userInfo);
    } catch (error) {
        errorResponce(error, res);
    }
};
export const logout = async (req: any, res: any) => {
    try {
        // const { email }
        // await userServices.deleteRefreshToken()
        res.status(200).send("Token is deleted");
    } catch (error) {
        errorResponce(error, res);
    }
};

// export const logout = async (req: any, res: any) => {
//     try {
//         const { rt } = req.header;
//         const userInfo = await userServices.clearToken(rt);
//         res.status(200).send(userInfo);
//     } catch (error) {
//         errorResponce(error, res);
//     }
// }

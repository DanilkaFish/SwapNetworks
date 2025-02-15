import { RefreshToken } from "@prisma/client";
import { hash } from "bcrypt";
import jwt from "jsonwebtoken";
// const jwt = require('jsonwebtoken');
require("dotenv").config();

const JWT_SECRET = process.env.JWT_SECRET;

export const SHA256 = (msg: string) => {
    return hash(msg, 10);
};

const SECRET_KEY = "alles";
// async function chech
export const generateRT = function () {
    return "toto";
};

export const generateJWT = function(payload: object) {
    const jwtsign = jwt.sign(payload, SECRET_KEY);
    console.log(jwtsign);
    return jwtsign;
};
export const jwtVerify = function(token: string) {
    return jwt.verify(token, SECRET_KEY);
}
// export const jwt = (header: any, payload: any) => {
//     return SHA256(String(header) + "." + String(payload) + "." + SECRET_KEY)
// };
